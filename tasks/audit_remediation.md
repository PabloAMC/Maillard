# Audit remediation tracker

Source: four-agent adversarial audit, 2026-08-26 (`scratch/gfi_call_brief.md`, untracked).
Branch: `audit-remediation` (off `s27-matrix-path-accuracy` @ 5fe8515 — strictly ahead of
`main`, so this branch supersedes both).

Legend: `[x]` done · `[~]` in progress · `[ ]` open · `[P]` needs Pablo (journal access or
science judgment) · `[D]` deferred by design (documented decision, not forgotten)

## Phase 1 — Data honesty (no behavior change)

- [x] 1.1 ProtocolPilot2026 files: relabel `origin` honestly (they are model output, not
      experiments), remove fabricated `measurement_date` / `uncertainty_pct`, state provenance
      in `notes`. Keep files so `hexanal_nonanal_calibration.py` lanes stay wired, but the
      lane semantics become "reproducibility", not "calibration_closed against experiment".
- [x] 1.2 Quarantine fabricated-source benchmarks to `data/benchmarks/quarantined/` (excluded
      from the panel by the loader's non-recursive glob): `cys_ribose_150C_Mottram1994`,
      `cys_glucose_150C_Farmer1999`, `thiamine_cys_ribose_100C_Hofmann1996`. Add README with
      the audit evidence and the nearest real candidate sources. Decision: quarantine, not
      delete — reversible, preserves the tolerance/values for forensics.
- [x] 1.3 Barrier provenance truth: `data/lit/arrhenius_params.yml` — flag the
      thiol_addition kJ/kcal collision (29 kJ/mol here vs 28.85 kcal/mol in
      `src/barrier_constants.py`) and the other YAML-vs-FAST divergences; relabel
      calibration-rationale comments in `FAST_BARRIERS` as fitted, not literature.
      Runtime values UNCHANGED (changing them silently re-calibrates the whole model).
- [x] 1.4 MFT odour threshold in `data/species/desirable_targets.yml` (0.0001 vs commonly
      cited ~0.005–0.007 µg/kg; README example table itself uses 0.007). Verify against
      literature, fix with provenance note.
- [x] 1.5 `data/lit/README.md:15` — `calibration_offsets.json` is an Optuna fit to two panel
      benchmarks (currently unwired), not MLP/xTB-vs-DFT offsets.
- [x] 1.6 Un-gitignore `results/validation/external_validation_report.*`; note in the report
      that the hold-out uses the wider `uncalibrated` prior tier than the in-panel headline.

## Phase 2 — Safe code fixes (test-backed, minimal behavior change)

- [x] 2.1 Guard: `calibrate_from_intake` (`src/matrix_experiment_intake.py`) must raise on
      `evidence_class == "external_validation_only"` so hold-out exclusion is enforced by
      code, not directory convention. + unit test.
- [x] 2.2 `src/benchmark_evaluator.py` water-activity key bug: fixed + tested. CORRECTION
      (2026-08-26 hygiene pass): the evaluator/reporting/assertions/registry/markdown module
      family is a DEAD-CODE ISLAND — nothing outside it imports it (it was even unimportable:
      BenchmarkThresholds imported from the wrong module, now repaired); the live path
      (benchmark_validation.py) always had the correct water_activity helper. So 2.2 fixed
      unreachable code; the divergence was never live. The island's twins are now at parity
      with sync comments; decision below on deleting vs reviving it.
- [x] 2.3 2-Pentylfuran SMILES in `benchmark_validation.py` / `benchmark_evaluator.py` is
      2-butylfuran (`CCCCC1=CC=CO1`); correct to `CCCCCC1=CC=CO1` (C9, MW 138.21) to match
      the curated registry. Behavior-affecting only inside the injection block (see 3.1).
- [x] 2.4 Heavy-atom-only MW (`sum(atom.GetMass())`) in `benchmark_validation.py:711` vs
      `Descriptors.MolWt` in `recommend.py` — use MolWt on both sides of the round-trip.

## Phase 3 — Science fixes (behavior-changing; batch, then regenerate once)

- [x] 3.1 mM-vs-M basis clash in the lipid-marker injection block
      (`benchmark_validation.py` ~L713: markers built in mol/L, dict consumed as mM →
      ~100× geomean suppression of Maillard targets on the network path). Minimal fix:
      convert to mM at injection. The deeper issue — injected markers are unreachable as
      outputs (`recommend.py` depth/exogenous filters), so the block currently contributes
      zero signal — documented; making them projectable is a design change for Pablo.
- [x] 3.x Physics guards (2026-08-26): protein_fraction>=0.999 sentinel in
      predict_headspace (defuses the 46x silent suppression landmine, clamp-and-warn
      mirroring recommend.py); van't Hoff temperature clamped to [273.15, 373.15] K
      (audit example: hexanal 453K gave Kaw 3.7; deliberately NO Kaw<=1 ceiling — real
      gases legitimately exceed 1). 22 new tests in test_thermo_basis_and_guards.py.
- [x] 3.y Joback entropy basis — AUDIT FINDING S8 OVERTURNED for balanced reactions:
      Joback formation-basis H and S carry the SAME frozen element offsets, which cancel
      identically in any atom-balanced reaction sum (proven algebraically + bit-identical
      numerics). Real exposure is only unbalanced steps: 197/198 SMIRKS steps are balanced
      (the one exception, thiamine Additive_Thermal_Degradation, sits below the gate), and
      thermodynamic gating resolves "off" for every panel benchmark anyway. Shipped: basis
      refactor with explicit delta_H/S/G + atom/charge-balance labelling (zero prediction
      change), and an opt-in NEUTRALIZE_UNBALANCED_THERMO_GATE flag (default False).
      [P] flipping that flag is more principled physics but removes an accidental
      mass-creation guard (test_thermo_gating relies on the element residual) — pair with
      an explicit upstream balance check if flipped.
- [x] 3.2 QRRHO free-rotor entropy (`src/thermo.py:75-89`): uses B_av as the moment of
      inertia and ignores `freq_cm1`; correct to Grimme 2012 reduced moment
      mu' = mu*B/(mu+B), mu = h/(8 pi^2 c nu). QM lane only. + unit test asserting the
      correction *caps* (not raises) low-mode entropy. Note: previously stored xTB/DFT
      barriers in results_db were computed with the wrong formula — flag for re-run.
- [x] 3.3 Regenerate validation artifacts after 1.2 + Phase 3 (panel shrinks; predictions
      move) and update tests that pin the old headline (e.g.
      `tests/scientific/test_lipid_oxidation_saturation.py` "37/48 guard").

## Phase 4 — Honest reporting

- [x] 4.1 README: replace 37/48 badge + calibration section with the split headline
      (literature rows vs internal reproducibility rows, reported separately, with median CI
      width) and the honest hold-out statement (0/8 genuine extrapolations; covered points
      re-score an in-panel anchor under the wider prior tier).
- [x] 4.2 Validation summary rendering: report literature and internal rows separately;
      print median CI width next to coverage; zero-width CIs → "not evaluable".
- [x] 4.3 VALIDATION_CONTRACT.md: add hold-out methodology section (prior tier difference).

## Deferred by design (documented, not forgotten)

- [D] Henry's-law observability gate (`recommend.py:296`) makes the partition layer inert and
      zeroes the `henry_kaw` MC channel. PREREQUISITE DONE (2026-08-26): the Kaw table was
      rebuilt against the full Sander 5.0.0 database (largest corrections: MFT 158x,
      acetone 136x, 2-pentylfuran 40x, acrolein 40x, DMTS 36x, FFT 24x; every entry carries
      previous_value + source + estimated flags; delta_H_sol refreshed; measured panel
      impact of all physics fixes combined: ONE unscored HMF row, 0/50 scored rows).
      Remaining: the gate flip itself, as a dedicated workstream with regeneration.
      [P] Two literature values deliberately NOT applied because they cross the 1e-8 gate
      and would flip observability classification: acrylamide 6e-8 (shipped 1e-9) and
      HMF 5e-8 (shipped 1e-10) — stored as literature_value_pending_threshold_decision;
      verified numerically inert today (acrylamide is produced by safety.py, not the
      gated loop). [P] FFT Kaw (2e-3) is the weakest entry — no measurement exists
      anywhere, ~5x uncertainty each way; matters the day the partition layer activates.
- [D] `sensory.py` air-vs-water ODT basis + squared retention — latent (prod passes
      `headspace=None`); fix alongside the Henry workstream.
- [D] "Limiting precursor" narrative vs geomean arithmetic in `projection.py` — decorative
      naming; revisit when the budget model is next touched.
- [x] Duplicate prediction function family (`benchmark_registry`/`benchmark_evaluator`/
      `benchmark_reporting`/`benchmark_assertions`/`benchmark_markdown`) — confirmed a
      closed dead-code island (stalled refactor; nothing imports it; one more latent bug
      left unfixed: benchmark_registry passes pH_profile= to ReactionConditions). Twins
      brought to parity + sync comments 2026-08-26. RESOLVED 2026-08-26 (owner-approved):
      island DELETED — re-verified that the only importers of the five modules were the
      modules themselves, so `benchmark_validation.py` / `presentation.py` remain the sole
      live prediction lane; the source-text assertions in
      `tests/unit/test_audit_remediation_2026_08.py` were updated to drop the dead module.

## Environment / owner-action items (2026-08-26 hygiene pass)

- [P] dvipng: needs `sudo tlmgr install dvipng` (TeX Live 2026basic is root-owned; no
      non-admin route). Then regenerate figures:
      `python scripts/generators/generate_validation_figures.py --docs-asset-dir docs/assets`
      and `generate_family_validation_figures.py --output-dir results/validation
      --docs-asset-dir docs/assets`. Stale PNGs: docs/assets/validation_overview.png (May 9),
      family_*.png (May 27) vs md/json refreshed Aug 26.
- [x] tests/qm NumPy/PySCF incompatibility: RESOLVED 2026-08-27 — numpy pinned to 2.3.5
      in the maillard env and environment.yml (>=2.0,<2.4 with rationale comment);
      tests/qm went 6 failed / 15 passed -> 0 failed / 21 passed / 3 skipped.
- [x] Docker parity: VERIFIED CLEAN 2026-08-27 — OrbStack started, existing container
      reused (no bootstrap), Docker-lane vs local-conda benchmark summary compared by
      direct snapshot diff (NOT git diff — the file is gitignored, so the originally
      written procedure was vacuous): markdown byte-identical; JSON differs only in
      mount-path strings and 4 last-ULP floats (~1e-16 relative). No scientific drift.
- [x] Root scripts: both tracked scratch scripts git mv'd into scripts/ingest/archive/
      2026-08-27 (local modifications preserved; archive README updated). Root is clean.

## Citation sweep results (2026-08-26, full ledger: results/validation/citation_verification_ledger.md)

- [x] Full anchor verification: 225 unique DOIs / 421 anchor sites. MATCH 88 (39%),
      METADATA-MISMATCH 47, TOPIC-MISMATCH 45, DEAD 45 — **61% of anchors defective**
      (the earlier 30-45% estimate was optimistic; the live-DOI-wrong-paper class is as
      large as the dead class). Damage concentrated in benchmark_intake_registry +
      slr_incorporation_matrix; executable benchmark JSONs mostly clean. 17 dead DOIs
      repaired in the registries (with doi_repair records); 9 arrhenius/prior flags added;
      duplicate-anchor collapse documented (single DOIs serving up to 4 "distinct" refs).
- [x] **Parker2012 quarantine review** — APPROVED AND EXECUTED 2026-08-26. Its DOI is dead
      with no identifiable source, the real Parker 2012 is a different system (finish-fried
      fries), and it carried the TIGHTEST tolerance in the collection (1.05/0.02 — the
      fitting tell); measured agreement was ratio 1.016, i.e. the model landed inside a
      ±5 % band on a value with no source. `git mv` to data/benchmarks/quarantined/ with
      source_metadata quarantine fields + README section. Its two audit-era test consumers
      were restructured, not re-retargeted (see Test fallout below).
- [P] **Hexanal ODT is unsourced**: 10.1021/jf00111a008 (cited for the off-note
      threshold) is actually a sweet-potato amino-acid paper. The 4.5 ppb value itself is
      consistent with common compilations — needs a real citation, likely not a value fix.
- [x] Cerny2008 benchmark DOI repaired 2026-08-26: source_doi 10.1021/jf800268t ->
      10.1021/jf801762c, with a source_metadata.doi_repair record {old, new, date, basis}.
      Values UNCHANGED — the file now also carries a prominent
      source_metadata.VALUES_NEED_RE_EXAMINATION warning that the paper's subject is
      5-hydroxy-3-mercapto-2-pentanone, not MFT, so the 2.47 ppb figure is not yet
      verified against the real paper and this benchmark must not constrain the sulfur
      barriers. [P] the values half still needs ACS access.
- [P] 28 dead DOIs unrepaired (4 are confabulation-signature non-DOIs); 6 benchmark-file
      repairs proposed-not-applied — see the ledger.

## Needs Pablo (journal access / science judgment)

- [x] Quarantined-benchmark verdicts — ALL THREE APPROVED AND EXECUTED 2026-08-26.
      (a) + (b) `git rm`'d (history preserves them); their README sections kept, marked
      "deleted 2026-08-26 after source-recovery confirmed no source exists". The lost
      chemistry regression is replaced by
      tests/scientific/test_pentose_hexose_sulfur_ordering.py, which asserts the ordinal
      literature constraint MFT(cys+ribose) >> MFT(cys+glucose) at matched conditions
      (10.1021/jf9705983) on synthetic formulations — observed 15.8x against a 3x floor —
      plus the SMIRKS-coverage assertions inherited from the deleted
      test_mottram_coverage.py / test_farmer_coverage.py. (c) rebuilt as
      data/benchmarks/thiamine_cys_glucose_120C_Bolton1994.json and IT FAILS HONESTLY:
      predicted 0.080 ppb vs the 13 ppb target = 163x under-prediction (MALE 2.21),
      robust at 35-745x across the assumed-value sweep. That is now the panel's real
      Family 03 signal. Registry entries repointed (see work log). Research record below
      retained verbatim:
      (a) Mottram1994: DELETE the ppb values — real Mottram 1994 (10.1021/bk-1994-0543.ch014)
      is norfuraneol+cysteine qualitative; 1995 paper is low-moisture+phospholipid,
      incompatible with the benchmark's aqueous a_w 0.95; no quantitative cys/ribose
      MFT+disulfide+furfural source exists (era reported FD factors / mol%). Fallback:
      ordinal pentose>>hexose MFT constraint from 10.1021/jf9705983.
      (b) Farmer1999: DELETE, no replacement — Farmer never published a cys/glucose model
      system (her corpus is cys+ribose); four chemistry tells mark the values invented;
      the real glucose/cys paper (10.1021/jf960456t) contradicts the benchmark.
      (c) Hofmann1996: REBUILDABLE as Bolton/Reineccius 1994 (10.1021/bk-1994-0543.ch022):
      thiamine+cys+GLUCOSE 120C/60min, MFT ~13 ppb, tolerance 3x, pH unpinned.
      Note: ACS Symp. Ser. 543 contains both real chapters — plausible common ancestor of
      both fabricated labels. calibrate_barriers.py is fitted to two nonexistent sources.
- [P] Re-anchor the WHOLE sulfur branch (widened from the thiol_addition units item,
      2026-08-26 forensics): commit 2ea7d12 ("Real tests based on experimental data",
      2026-03-15) set thiol_addition/thiol_addition_hexose/thiol_oxidation/
      aminoketone_condensation/strecker_degradation in one Optuna-style fitting loop
      whose explicit targets included the now-quarantined Mottram1994 + Farmer1999
      (documented verbatim in that commit's todo/lessons diffs). Post-quarantine
      constraint status: thiol_addition is the ONLY one still literature-constrained —
      by Hofmann1998 alone, to a 0.75 kcal window [28.10, 28.85], and the shipped 28.85
      sits at the boundary (Hofmann-optimal is ~28.60, MALE 0.082 -> 0.051). The others
      have NO surviving literature constraint (their observables are quarantined or
      synthetic) and are free parameters. The xylose/thiamine benchmarks do NOT
      independently confirm 28.85 — their per-lane upstream_observability_factor was
      fit with the barrier frozen. Also: Hofmann's own tolerance (1.45/0.09) is the
      3rd-tightest in the collection — same fitting-tell pattern. cysteine_thermolysis
      is inert (0/16 benchmarks move on +-3 kcal): its 206-vs-125 kJ/mol dispute is
      documentation debt only.
- [P] cysteine_thermolysis anchor: RESEARCHED 2026-08-26 — best real anchor is
      Zheng & Ho 1994 (10.1021/bk-1994-0564.ch011): aqueous cysteine H2S release,
      Ea 123-135 kJ/mol across pH 3-9 (mean 130). FAST_BARRIERS' 30 kcal/mol is thus
      ALREADY right; the YAML's 206 kJ/mol is the dry-pyrolysis regime, and its A/Ea
      pair is internally broken (k(100C)~5e-22/s — inert on cooking timescales; a
      consistent A at Ea~130 is ~1e14-1e15/s). Runtime inert anyway (forensics: 0/16).
      Decision: rewrite the YAML entry to Zheng & Ho.
- [P] Martins Schiff Ea: RESEARCHED 2026-08-26 from the actual Wageningen thesis
      (edepot.wur.nl/121418, Table 5.2) — glucose+glycine condensation step reports
      Ea = 97 +/- 3 kJ/mol at pH 6.8; 57 kJ/mol appears NOWHERE in the thesis, and the
      A_value comment is not literal either. Caveat: 97 is a condensation+Amadori
      composite — raising this entry without re-examining `amadori` double-counts.
      Decision: correct to 97 with the composite caveat, or re-derive the split.
- [P] Whether to make injected lipid markers projectable observables (3.1 design question).
- [x] Bulk backlog DOI resolution (2026-08-26): coverage 53/239 -> 138/239 (65 high +
      20 medium confidence, hand-reviewed; 32 self-referential pointers machine-tagged;
      69 left for lack of signal; 2 contradicting identifiers deliberately NOT written).
      [P] the 69 no-signal items + the 2 flagged contradictions need human eyes.
- [x] Sigma derivation (2026-08-26): leave-lane-out residuals give RMS ln-sigma 2.86,
      90% CI [1.98, 5.48], n=6 — shipped 2.0 is at the LOWER edge (conservative, not
      hold-out-widened). Artifact: results/validation/matrix_sigma_residual_derivation.md.
      [P] decide whether to raise sigma to ~2.9 (would widen uncalibrated-tier CIs).
- [P] Furfural temperature turnover: RESOLVED mechanistically (2026-08-26) — it is a
      projection-layer ARTIFACT, not chemistry. Allocation set is constant (5 species,
      110-250C); the cause is a saturating budget (severity sigmoid centered 110C/width
      18 caps the pool at 1.108x its 150C value) divided by a still-growing softmax
      denominator, plus a double-applied Boltzmann on path_span (recommend.py:743 vs
      745-748, effective T/2.19). Furfural's activity is pinned at 1.0 — its output is
      pure normalization. Real Arrhenius drive 150->170C is 4.5x vs budget growth 1.07x.
      Decision needed: retune _projection_temperature_factor (projection.py:52-64) to
      track Arrhenius, and de-duplicate the Boltzmann — both change every free-precursor
      prediction (recalibration event).
      RESOLVED AND SHIPPED 2026-08-27 (owner-approved, final reconciliation wave).
      (a) The severity sigmoid is DEMOTED, not deleted: it survives as the bounded
      process-state index consumed by melanoidin trapping / depth bias / direct-sulfur
      bonus (heuristics that should saturate) and no longer scales the budget.
      (b) The budget's thermal dependence is now a first-order conversion extent under an
      apparent Arrhenius rate: drive = (t/tau_ref) * exp(-(Ea/R)(1/T - 1/T_ref)),
      yield = baseline + (ceiling - baseline) * (1 - exp(-drive)), with ceiling pinned at
      1.0 (mass conservation on the limiting pool), Ea pinned at 120 kJ/mol (literature:
      the `enolisation` entry in arrhenius_params.yml, inside the 97-160 envelope spanned
      by the other volatile-forming steps) and T_ref 423.15 K.
      (c) The double Boltzmann is gone: `span_activity` (a second exp(-dspan/(0.65 RT)))
      and the ^0.65 softening on the flux were BOTH removed; the allocation now uses the
      relative pathway flux directly, i.e. exp(-span/RT) exactly once at physical T,
      which is also the correct steady-state partition k_i/sum(k_j).
      (d) Two free constants refit against LITERATURE ROWS ONLY (synthetic + hold-out
      excluded by assertion, not convention): baseline_volatile_yield_fraction 1.0e-6
      (unidentified — flat to 1e-9, incumbent kept) and reference_conversion_time_min
      1.2589e4 (10**4.10). Reproducible: scripts/generators/refit_projection_constants.py,
      record in results/validation/projection_constant_refit.{json,md}.
      MEASURED COST, reported not tuned away: literature-row MALE 0.147 -> 0.261 dex over
      all 24 rows, and 0.309 -> 0.558 dex (0.119 -> 0.390 excluding the known Bolton
      failure) over the 11 rows the projection layer actually moves. Essentially all of
      the gap is the Boltzmann de-duplication — the budget retune is a single global scale
      that the tau_ref fit absorbs. That gap is the measure of how much of the panel's
      previous agreement on multi-analyte free-precursor benchmarks was being carried by
      a selectivity term evaluated at less than half the physical temperature.
      Furfural is now MONOTONE in temperature (62 / 352 / 1688 / 7154 / 25770 ppb at
      110 / 130 / 150 / 170 / 190 C); the peak is gone for the right reason. Pinned by
      tests/scientific/test_time_propagation.py::
      test_furfural_output_is_monotonic_in_temperature_across_the_full_range and
      tests/unit/test_budget_projection.py::
      test_projection_budget_temperature_dependence_tracks_arrhenius_not_a_sigmoid.
      [P] NOT refit, flagged for a dedicated workstream: recommend.py
      `depth_bias_strength` (0.85 offset) and `direct_sulfur_bonus` (0.8 coefficient) are
      unconstrained legacy fits from the quarantined-target era. direct_sulfur_bonus in
      particular is the exact knob that would absorb the MFT residuals this panel exists
      to expose, so tuning it now would launder the finding above.

## Work log

- 2026-08-26: audit completed (4 agents); branch `audit-remediation` created off s27;
  tracker created.
- 2026-08-26 (session 2): Phases 1-3 executed. 1.1: ProtocolPilot intake YAMLs +
  JSONs relabeled synthetic_diagnostic/diagnostic_only, fabricated dates and
  uncertainties removed, lanes artifact carries an evidence disclaimer; new
  generator scripts/generators/refresh_internal_reproducibility_snapshots.py
  re-freezes Internal2026 + ProtocolPilot snapshots after intentional model
  changes (replaces the previous hand-frozen commits). 2.x/3.1/3.2: all landed
  with tests in tests/unit/test_audit_remediation_2026_08.py (5 passing); the
  injection fixes changed matrix_precursor_augmented predictions, so snapshots
  and their ranking contracts were re-frozen (expected drift, documented in the
  snapshot notes). Also fixed a dict(None) crash in primary_benchmark_campaign
  when no promotion target qualifies. 4.2: prediction_uncertainty payload +
  markdown now report a signal-origin coverage split (external_literature vs
  internal_synthetic, zero-width rows -> not_evaluable, median CI widths). 4.3:
  VALIDATION_CONTRACT gained the provenance-audit banner, hold-out methodology
  (S3E) and coverage-split reading rules (S3F); stale strict-ready claims for
  quarantined benchmarks withdrawn. 1.2/3.3 delegated to a subagent (quarantine
  moves, ~17 test files, artifact regeneration). Known environmental failure:
  test_cli_scripts.py::test_run_campaign_named_comparison_mode needs dvipng
  (LaTeX plotting), pre-existing on this machine.
- 2026-08-26 (session 2, cont.): quarantine executed by subagent — 3 benchmarks moved
  to data/benchmarks/quarantined/ with forensic README + registry annotations; ~17 test
  files updated; artifacts regenerated locally (prediction_uncertainty n=200 ~6 min;
  external_validation_report now tracked; validation_overview md/json refreshed; PNGs
  STALE pending Docker/dvipng). New honest numbers: panel 16 benchmarks / 41 matched rows;
  aggregate 35/41; split = external_literature 9/10 evaluable in CI (median width 3.13 dex,
  5 not-evaluable) vs internal_synthetic 18/18 (4.25 dex, 8 not-evaluable); hold-out
  nominal 3/8 = 3/3 re-scoring bundles + 0/5 genuine extrapolations. strict-ready 8/16.
  Post-quarantine semantics: new external_data_status "synthetic_diagnostic_only" (maps to
  reference-grade support, honest promotion blocker text); campaign resolver accepts
  synthetic origin; example intake contract ranks refreshed to current model. README 4.1
  rewritten (yellow calibration badge, split table, honest hold-out + provenance notes).
  Full-suite verification run pending; known environmental failure: dvipng CLI test.
- 2026-08-26 (session 2, round 2 close): five parallel streams landed — defensibility
  math (sigma derivation, furfural artifact diagnosis, sulfur-tuning forensics), full
  citation sweep (61% anchors defective; 17 DOI repairs; tracked ledger), source recovery
  (Mottram/Farmer = delete verdicts, Hofmann1996 rebuildable as Bolton 1994; Schiff Ea
  97 kJ/mol from the thesis; Zheng & Ho vindicates FAST cysteine value), physics fixes
  (Kaw table vs Sander 5.0.0; guards; Joback finding overturned for balanced reactions),
  hygiene (dead-code island found + twins parity; scripts archived; qm-lane numpy/pyscf
  breakage identified; dvipng + Docker parity = owner actions). Central regeneration
  re-run: aggregate 35/41 and the 9/10 / 18/18 split reproduced exactly — physics fixes
  moved zero scored rows, as measured. Evidence artifacts (sigma derivation, citation
  ledger) git-tracked. All work remains uncommitted pending owner instruction.
- 2026-08-26 (round 3, owner approved all decisions): main-session value edits applied —
  thiol_addition 28.85 -> 28.60 (Hofmann-only optimum); thiol_oxidation/hexose/
  aminoketone/strecker relabeled UNCONSTRAINED LEGACY FIT (values kept for continuity);
  uncalibrated matrix_headspace sigma 2.0 -> 2.86 (residual-derived point estimate,
  ~+/-110x at 90%); acrylamide Kaw 1e-9 -> 6e-8 and HMF 1e-10 -> 5e-8 activated
  (previous values retained in-file); arrhenius YAML schiff Ea 57 -> 97 kJ/mol (thesis
  Table 5.2, composite caveat noted vs amadori) and cysteine_thermolysis re-fitted as a
  pair to Zheng & Ho (Ea 130.4 kJ/mol, A 1e14/s); README + VALIDATION_CONTRACT sigma
  texts updated. Barrier-adjacent unit tests green; scale-pinned scientific tests will
  drift until the final reconciliation wave. In flight: Parker2012 quarantine + Bolton
  rebuild; island deletion + numpy pin + parity; two reference-repair partitions.
  Pending: projection budget retune + Boltzmann de-dup + central regeneration + full
  test reconciliation (final wave).
- 2026-08-26 (session 2, panel surgery — owner-approved decisions executed): four decisions
  from "Citation sweep results" + "Needs Pablo" landed as a single panel-surgery pass.
  (1) `acrylamide_asparagine_glucose_Parker2012` QUARANTINED (git mv + source_metadata
  quarantine fields + README section): dead DOI, no identifiable source, tightest contract
  in the collection (1.05/0.02) with an observed ratio of 1.016 — the fitting tell.
  (2) `cys_ribose_150C_Mottram1994` + `cys_glucose_150C_Farmer1999` DELETED (git rm) after
  source recovery; README sections kept as the record. Replacement regression:
  tests/scientific/test_pentose_hexose_sulfur_ordering.py (ordinal pentose >> hexose MFT,
  10.1021/jf9705983; observed 15.8x vs a 3x floor) which also absorbs the SMIRKS-coverage
  assertions from the now-deleted test_mottram_coverage.py / test_farmer_coverage.py.
  (3) `thiamine_cys_ribose_100C_Hofmann1996` DELETED and REBUILT from the verified Bolton
  1994 chapter as `data/benchmarks/thiamine_cys_glucose_120C_Bolton1994.json`
  (10.1021/bk-1994-0543.ch022; 120 C/60 min; MFT 13 ppb; tolerance 3.0x/0.48; pH 5.0 and the
  molarities explicitly marked assumed_not_from_source). It evaluates supported, coverage
  1.0, predicted 0.080 ppb vs 13 ppb = 163x UNDER-prediction (MALE 2.21) — a genuine,
  robust Family 03 failure (35-745x across the assumption sweep), replacing a fabricated
  row that "passed". Registry entries `hofmann_schieberle_grosch_1996` in
  benchmark_intake_registry.json and slr_incorporation_matrix.json repointed to the Bolton
  source/benchmark with `source_repair` records; the ids were kept as stable handles
  (computational_priors.json, deep_research_backlog.json and a test key off them) and both
  entries now say so explicitly. The old anchors "conversion 0.021%" / "Cys 2.3x MFT" were
  WITHDRAWN in the SLR matrix as unsourced. (4) Cerny2008 DOI repaired with the
  values-need-re-examination warning (see the item above).
  [P] LEFT ALONE DELIBERATELY: `data/lit/computational_priors.json` still carries an entry
  id `hofmann_schieberle_grosch_1996` with `confidence_tier: "high"` and the two withdrawn
  numbers (thiamine_to_MFT_conversion_pct 0.021, cysteine_augmentation_MFT_fold 2.3), and it
  is consumed by src/literature_runtime.py. Those figures came from the fabricated source.
  Correcting them is a runtime-behaviour change and was out of scope for this pass — it
  needs an owner decision alongside the sulfur re-anchoring item.
  DONE 2026-08-27 (final reconciliation wave, Part 2b): downgraded honestly.
  confidence_tier "high" -> "low" (the lowest tier in use in that file);
  provenance_tier -> "unsourced_withdrawn"; uncertainty_posture -> "withdrawn_do_not_use";
  new source_status "no_verifiable_source"; new source_repair record and a notes field
  pointing at the replacement lineage (Bolton & Reineccius 1994,
  10.1021/bk-1994-0543.ch022 -> data/benchmarks/thiamine_cys_glucose_120C_Bolton1994.json).
  The two numbers are RETAINED in-file (not deleted) with a `withdrawn_values` list, so the
  forensic record survives and the id keeps working as a join key.
  CONSUMER AUDIT (the question this item asked): the two numeric fields
  `thiamine_to_MFT_conversion_pct` and `cysteine_augmentation_MFT_fold` have NO live
  consumer — the only references anywhere in the repo are this file and an archived ingest
  script (scripts/ingest/archive/ingest_batch_alternative_matrix.py). The entry does reach
  src/literature_runtime.py, but only as one row in
  `query_family_runtime_priors(runtime_lane="thiamine_fragmentation")` listings; it is NOT
  the lane default (that is `cerny_2007_thiamine_split_v1` via
  `_RUNTIME_LANE_PRIOR_REGISTRY`), and `get_family_runtime_prior` therefore never selects
  it. Sensitivity to the two withdrawn numbers is exactly ZERO: no prediction, benchmark
  row or artifact number moves. The downgrade changes only what the entry advertises about
  itself in prior listings.
  Panel arithmetic for the next regeneration: size UNCHANGED at 16 (-Parker2012,
  +Bolton1994), row count unchanged (both files carry exactly one analyte), but the swap
  costs one PASS (Parker acrylamide, ratio 1.016, strict-ready) and adds one FAIL (Bolton
  MFT, ratio 163). Do NOT hand-adjust the old 35/41 + 8/16 headlines: other in-flight work
  on this branch moves scored rows too (thiol_addition re-centred 28.85 -> 28.60; the
  Sander 5.0.0 Kaw refresh has HMF/acrylamide crossing the 1e-8 observability gate, so the
  panel now has ZERO low_headspace rows). Live measurement at the end of this pass, for
  reference only: 16 benchmarks, 50 comparison rows, 43 within 2x, 6 strict-ready. The
  headline must be recomputed by the central regeneration, not predicted. No validation
  artifacts were regenerated in this pass, by instruction.
  [P] NEW REPORTING BUG SURFACED (pre-existing, newly consequential):
  scripts/generators/generate_validation_figures.py computes `worst_quantitative_ratio`
  over ALL quantitative benchmarks but `worst_quantitative_point` only over the
  reference_volatiles subset (L162: `max(reference_quantitative_benchmarks, ...)`). With
  Bolton in the panel the two now disagree loudly -- the artifact will report "worst
  quantitative benchmark ratio: 162.754x" next to a "worst point" of Cerny at 1.641x. Until
  now the two happened to coincide, so the naming mismatch was invisible. Left unfixed here
  because it changes report output; fix it in the regeneration workstream (either rename the
  key to worst_reference_quantitative_point or widen its population).
  FIXED 2026-08-27 (final reconciliation wave, Part 2a): population WIDENED.
  `worst_quantitative_benchmark` is now selected over all quantitative benchmarks, so it
  pairs with `worst_quantitative_ratio` (both now report Bolton at 173.315x), and a new
  `worst_quantitative_population` key states the population explicitly. The reference-only
  view it used to provide is preserved under `reference_worst_quantitative_ratio` /
  `reference_worst_quantitative_point` (Cerny, 3.419x), and both are rendered as separate,
  labelled lines in validation_overview.md.
  `tests/scientific/test_validation_figures.py` had PINNED the buggy pairing; it now
  asserts the consistency property (point.max_ratio == worst_quantitative_ratio) instead
  of a hard-coded benchmark id, so the pair cannot silently diverge again.
- 2026-08-27 (reference partition B: SLR/priors/backlog/payload files): ALL defects
  dispositioned — 103 ledger-flagged DOIs (55 repaired/corrected, 19 no_verifiable_source,
  3 identifier_unavailable, rest confirmed-correct or escalated+fixed), 113 deep-research
  refs (67 repaired, 43 no_verifiable_source, 3 unresolvable), backlog now has 0
  undispositioned items. 214 DOI fields in the partition, 0 dead. 94 shared_anchor
  annotations over 25 duplicate DOIs. Zero numeric leaves moved (verified vs HEAD).
  Highlights: several stored values CONFIRMED verbatim in the true sources (Sun 2015
  CML/CEL Ea pair; Yu 2020 EGCG kinetics); several exposed as unanchored or wrong
  (Ohsu kokumi EC50s three orders of magnitude off + out-of-range sensory units;
  mottram_2001 MFT-quench anchor has the WRONG SIGN vs real carnosine literature
  (+61% stored, literature reports a decrease); siripitakpong FFT retention anchor is
  actually vanillin data). 91 citation_caveat records mark right-paper-wrong-number
  cases. [P] NEW science calls: the 43 unanchored-value refs (worst: R5P donor
  hierarchy, HDMF anchor, beta-CD K_assoc table, MFT-quench sign), the vanillin-vs-FFT
  anchor, 5 fabricated-author paper_ids kept as join keys (id_provenance_note).
  Hexanal ODT: true source is a no-DOI Grosch 1982 book chapter, value unconfirmed —
  still open in partition A's hands.
- 2026-08-27 (readiness phase, owner approved "all of the above"): launched in parallel —
  (1) adversarial hunt over un-audited modules (safety.py first, then optimizer/extrusion/
  sensory-live/kokumi/pathway/results_db/input layers); (2) SMIRKS chemistry soundness
  review (mechanism-level, per template/family, cys+ribose and cys+glucose networks);
  (3) front-door smoke: README/QUICKSTART commands run verbatim + report-structure honesty
  + determinism (same-seed double run); (4) CI authoring (.github/workflows/ci.yml + an
  offline citation-contamination gate + a hold-out-exclusion guard script). Sequenced
  next: projection retune + central regeneration + full test reconciliation (after the
  Parker/Bolton agent lands), then stale-docs sweep + README sample refresh + AUDIT.md,
  then structured commits, then a zero-context cold-start red team against the settled
  tree as the final verification.
- 2026-08-27 (reference partition A: intake registry + benchmarks + species): 125
  ledger-flagged ids dispositioned — 98 repaired (56/56 replacement DOIs CrossRef-verified
  before writing), 24 no_verifiable_source, 3 left with documented reasons. Standouts:
  hexanal ODT re-anchored to the Belitz/Grosch/Schieberle compilation (both partitions'
  divergent calls harmonized by main session in flavor_reference_payloads.json, caveat
  retained: 4.5 ppb is compilation-level, not primary-table-verified); a DOI collision
  parasitising the Bi 2020 hold-out anchor un-collided; a swapped DOI pair fixed; four
  PMC/PMID entries recovered via embedded ids; three fabricated identifiers exposed
  (repo ids pasted into doi fields). Zero model-consumed numbers changed (verified).
  [P] ramirez_jimenez_2000 carries matrix_family "mild_legume_extrudate" on a toasted-
  BREAD study — the field is consumed by routing code, annotated not rewritten; owner
  call. Combined across both partitions: every flagged anchor in the repo is now either
  repaired with a verified source or explicitly labeled no_verifiable_source.
- 2026-08-27 (panel surgery landed): Parker2012 quarantined (predicted-vs-"measured"
  ratio 1.016 against a 1.05 contract — the fitting tell quantified); Mottram1994 +
  Farmer1999 + Hofmann1996 deleted (git history preserves; README sections retained);
  Bolton1994 thiamine benchmark created (verified source, MFT 13 ppb, 3x tolerance) and
  the model HONESTLY FAILS it at 163x under-prediction (sensitivity-swept 35-745x across
  assumption ranges; inverted test pins the failure). New ordinal pentose>>hexose MFT
  test passes at 15.8x vs 3x floor. Panel stays 16 benchmarks. NEW items for the final
  wave: (a) generate_validation_figures.py computes worst_quantitative_ratio over all
  benchmarks but worst_quantitative_point over reference_volatiles only — inconsistent
  pair with Bolton in panel; (b) computational_priors.json entry
  hofmann_schieberle_grosch_1996 still carries confidence_tier high + the withdrawn
  fabricated numbers (0.021%, 2.3x) and is consumed by literature_runtime.py.
- 2026-08-27 (SMIRKS chemistry soundness review — MAJOR): mechanism-level review of all
  ~35 emitted reaction families. Core sugar trunk (Schiff/Amadori/3-deoxyosone/furfural/
  HMF/Strecker map/beta-elim/CML-CEL) is SOUND. But: (1) MFT — the flagship — is made by
  a fabricated one-step 3-deoxyosone+2H2S shortcut; the accepted route (2,3-enolisation
  -> 1-deoxyosone -> norfuraneol -> +H2S/reductone -> MFT, van den Ouweland & Peer 1975)
  is entirely absent, and ribose-beats-glucose for MFT rests on a 1.05 kcal difference
  between two unconstrained fitted barriers; (2) fictitious H2 economy: free H2 in five
  families incl. a closed loop (Strecker+deamination) that manufactures H2 from water and
  feeds FFT production; (3) atomic [S] as a balance token, with Radical_Crosstalk
  consuming MFT to mop it up; (4) template-layer bis(2-methyl-3-furyl) disulfide is the
  WRONG REGIOISOMER (curated layer is right — two layers, two molecules, one label);
  (5) DMHF/HEMF both wrong isomers; pentose->DMHF unbalanced by one carbon; norfuraneol
  missing; (6) lipid radical chain NON-FUNCTIONAL: radical_propagation_o2 matches any
  sp2 carbon (61/103 steps are fabricated peroxidations), radical flags lost in
  _apply_smirks_rule ordering, beta-scission pattern requires sp3 beta-C so HEXANAL IS
  STRUCTURALLY UNREACHABLE in the network (masked in prod by the fixed branching ratio);
  (7) EIGHT families silently disabled by the 45 kcal/mol DEFAULT_BARRIER fallthrough
  (acrylamide, CML/CEL, furanone, thiamine, GSH lanes = kinetically off, half-life
  ~39,000 years); (8) README's "16 families of reaction chemistry" conflates SLR lanes
  with reaction generation — ~6/16 have genuine templates, 5 have none at all; the
  family-13 DFT prior is misapplied to furfural+H2S; curated pathway C is not
  self-contained and curated barrier wiring keys on family names the engine never emits.
  CLEAN LIST recorded (12 sound families). PLANNED WAVE G (after the reconciliation
  wave lands, to avoid refit-against-moving-target): mechanical fixes — radical-explicit
  SMIRKS, _fix_radicals ordering, correct furanone/disulfide SMILES, glyoxal-for-
  glycolaldehyde in thiazole condensation, explicit FAST_BARRIERS for the 8 fallthrough
  families (labeled estimated — note this ACTIVATES safety lanes), README family-table
  honesty split (SLR lanes vs implemented chemistry) — plus the deep fix: implement the
  real 1-deoxyosone/norfuraneol MFT route with honestly-tiered barrier priors, refit the
  constrained knob on Hofmann1998 only, regenerate, reconcile.
- 2026-08-27 (CI landed): .github/workflows/ci.yml added (three tiers: 1-min dependency-
  free gates job -> tests-unit -> tests-full; LaTeX installed in CI rather than
  deselected; qm behind workflow_dispatch; env built the Dockerfile's proven way with
  CPU-torch pre-install; numpy pin asserted). scripts/ci/citation_gate.py (confabulation-
  signature + syntax + repair-record + status-coherence checks, ratcheting baseline,
  negative-tested, zero false positives over 236 DOIs) and scripts/ci/holdout_guard.py
  (caught 5/5 injected regressions) both PASS. Gate findings acted on: 9 stale
  no_verifiable_source flags cleared in the backlog (repair pass had left pre-repair
  flags) and the gate baseline ratcheted down to 4 waivers — the 4 real DOI-less
  sources (2 theses, 1 patent, 1 report) sitting in doi fields. [WAVE G] give those a
  typed identifier/identifier_scheme pair (note: Liu hold-out bundle's source_doi
  truthiness feeds _matrix_external_data_status — handle with care). Known: first real
  Actions run may need tweaks (defaults-channel ToS, solve time); gates job is
  dependency-free and green immediately.
- 2026-08-27 (un-audited module hunt — MAJOR, defines Wave G scope with the SMIRKS
  review): TIER 1 SAFETY: (1.1) CML/CEL/furosine predictions consumed as ppb in the
  benchmark lane but as mg/kg / mg/100g-protein in literature_runtime — the "validated"
  safety rows (ratios 1.005-1.204) are an artifact of a unit collision (32 mg/kg written
  as conc_ppb 32.0 = 1000x); model CML is ~1200x below the literature band it claims to
  match; (1.2) molar_ratios have NO declared unit anywhere — acrylamide ppb varies 3300x
  across the three authoring conventions found in-repo; (1.3/1.4) aggregate risk score
  is log-compressed (1e9 span in acrylamide moves it 0.9), never in its documented [0,1],
  and the 100-ppb threshold + SafetyResult.flagged are dead code (nothing regulatory is
  consumed; the EU "meat analogue" benchmark cited does not exist in Reg. 2017/2158);
  (1.5) cys/gly mitigation is presence-only (trace = 65% discount); (1.6) runtime safety
  lane hard-pins pH=6.0 (86x low at pH 8); (1.7/1.8) furosine/LAL flagged from ingredient
  alone at 25C, 0.01C step flips flags at 121C; (1.14) "asn" substring parses asparaginASE
  as 20 mM asparagine. TIER 2 OPTIMIZER: (2.1) CLI's reported best formulation discards
  every optimized concentration (Optuna param names never match precursor labels -> flat
  1:1:1) and trial pH leaks into L-Phenylalanine's concentration; prints mM as "M";
  (2.2) risk_aversion is inert (safety constant per study); (2.3) grid ranking and BO
  optimize DIFFERENT objectives and disagree on the shipped grid; (2.4) wheat_gluten CLI
  choice always crashes; (2.7) explicitly disabling yeast_fermentation turns it ON.
  TIER 3: (3.1) every SME quantity saturated across the real 300-800 kJ/kg window;
  (3.2) the wet-lab DoE proposes two arms the model predicts byte-identically; (3.3)
  pathway ranker sorts on span only — pH/T/aw rates computed then ignored; (3.7)
  results_db SILENTLY OVERRIDES the frozen FAST barrier (Thiol_Addition_H2: DB 15.0 vs
  curated 28.6 = 1e7x rate) and the answer depends on the process CWD (relative path);
  (3.8) method preference inverted (xTB beats DFT; DFT rows get the larger sigma into
  published CIs); (3.9) 11 pre-QRRHO-fix DFT rows ship in the tracked DB, dormant only
  because of a second bug — DO NOT fix the reaction_meta key without purging/versioning
  first. TIER 4 SENSORY: (4.1) beany->meaty masking inert on 8/12 meaty compounds (list-
  vs-substring semantics); (4.2) 2,5-DMP synergizes with itself (+30%); (4.3) fabricated
  0.1 ppb ODT default assigned to all 8 toxic markers (latent); (4.5) QDA "0-10 scale"
  is self-normalized (top note always 10.0). EXTENSIVE CLEAN BILL recorded per module
  (acrylamide kinetics math exact; ranker span arithmetic fuzz-verified 2000/2000;
  precursor_resolver exemplary; determinism bit-exact). Six-fix priority list captured.
- 2026-08-27 (front-door smoke + determinism): SEVERE doc-vs-code drift found.
  (1) The README headline command silently DISCARDS the user's formulation: --target
  routes to inverse design over the fixed 15-entry grid; --sugars/--amino-acids/--ratios/
  --time-minutes are never read (proven: changing sugars moves 2/6202 leaves, ULP noise),
  yet the report's SS1 echoes the user's args while the science comes from "Cysteine
  Enrichment (Basic)". Forward mode honours the formulation but is undocumented.
  (2) README Step 2 (optimize_formulation) ALWAYS crashes at the end: best_trial.params
  includes categorical strings coerced via float() in build_family_upstream_contract —
  same root as module-hunt 2.1. (3) Quickstart comparison references grid names that
  don't exist ("Baseline"). (4) Ingest's documented flags omit two required fields;
  preview mode writes 4 files into results/validation despite no --confirm; README's
  "updates the benchmark database and regenerates artifacts" claim is false (writes one
  YAML). (5) README sample-output table's Confidence vocabulary (bounded_calibration/
  transferred_literature/surrogate_family) is FICTIONAL at compound level — the real
  column emits exploratory/low/medium/high, and THREE incompatible tier vocabularies
  coexist in one report.json (tier, confidence_tier, and the report's own SS6 glossary
  which matches none of its own columns). (6) Determinism: scientifically deterministic;
  bit-identical ONLY under PYTHONHASHSEED=0 (hash-order-dependent summation, ~1 ULP);
  no --seed flag; undocumented. (7) generate_validated_envelope_report.py is the only
  generator without --output-dir. (8) data/campaigns/ referenced by QUICKSTART does not
  exist. (9) GLOSSARY: 7/13 defined terms never appear in output; ~15 emitted terms
  undefined. (10) reporting.py:1758 emits a malformed 10-field row in an 8-column table.
  All -> WAVE G scope (docs/UX/tier-vocabulary unification stream).

- 2026-08-27 (FINAL RECONCILIATION WAVE — projection retune, Part 2 fixes, central
  regeneration, test + text reconciliation). Five parts, all landed.
  PART 1 (projection budget retune + Boltzmann de-dup): see the fully-expanded record
  under "Furfural temperature turnover" above. Headline: budget thermal dependence is now
  an apparent-Arrhenius first-order conversion extent (Ea 120 kJ/mol, literature-pinned;
  ceiling 1.0 by mass conservation; T_ref 423.15 K); Boltzmann selectivity applied ONCE at
  physical T (was T/2.19); two free constants refit against literature rows only via the
  new, reproducible scripts/generators/refit_projection_constants.py. Literature-row MALE
  0.147 -> 0.261 dex — WORSE, and reported as the finding rather than tuned away. Furfural
  monotone 110-190 C. Synthetic reproducibility snapshots refreshed (generator tag v4).
  PART 2: (a) generate_validation_figures.py worst-quantitative population mismatch FIXED
  by widening; reference-only view preserved under explicit new keys. (b)
  computational_priors.json `hofmann_schieberle_grosch_1996` downgraded to
  confidence_tier "low" + source_status "no_verifiable_source" + replacement lineage note;
  consumer audit found the two withdrawn numbers have ZERO live consumers and the entry is
  not the thiamine lane default, so sensitivity is exactly zero.
  PART 3 (central regeneration, tracked artifacts): benchmark_summary, benchmark_index,
  prediction_uncertainty (n=200, ~6 min), external_validation_report,
  experiment_value_ranking, loo_leverage, gap_heatmap + experiment_brief_cards,
  experiment_requests (+ stale rank-suffixed files pruned), validation_overview md/json and
  family_validation_overview md/json via non-figure drivers (PNGs still blocked by the
  missing dvipng), matrix_sigma_residual_derivation, projection_constant_refit.
  Two artifact-honesty bugs found and fixed during regeneration:
  derive_matrix_sigma_from_residuals.py hard-coded `shipped = 2.0` and
  src/external_validation.py hard-coded "ln-sigma 2.0, ~±27x" in the report's methodology
  disclosure — both silently wrong since the sigma was raised to 2.86 on 2026-08-26. Both
  now read the live constant, so the artifacts cannot drift from the runtime again.
  NEW HEADLINE NUMBERS: panel 16 benchmarks / 41 matched rows / 13 in the MC panel;
  aggregate CI coverage 33/41 (80.5%); split = external_literature 7/11 (63.6%, median CI
  width 1.97 dex, 4 not-evaluable) vs internal_synthetic 18/18 (100%, 2.66 dex, 8
  not-evaluable); hold-out 5/8 (62.5%) = 3/3 re-scoring + 2/5 genuine extrapolation, median
  fold error 32.79x (was 33.84x — the coverage gain is the wider sigma, NOT better
  predictions); strict-ready 4/16 (was 6/16: spi_hvp and wheat_gluten lost it to the
  retune); benchmarks without blocking gaps 11/16; out-of-CI cells 8/41.
  PART 4 (test reconciliation): full suite green except the known dvipng failure. Every
  changed expectation carries a dated comment naming its cause; nothing was silently
  relaxed, and no benchmark JSON contract tolerance was touched, so all the honest failures
  stay visible in the regenerated artifacts.
  PART 5 (text coherence): README calibration badge + split table + hold-out note +
  parity/experiment-priority rows updated to the regenerated values, with an added
  paragraph stating plainly that the calibration number FELL because of the retune and why;
  VALIDATION_CONTRACT §2 strict-ready delta and §3E current hold-out numbers added;
  CURATING_LIPID_OXIDATION_ANCHOR hold-out line updated with the wider-prior caveat.
  [P] CARRIED FORWARD, not addressed here (out of this wave's scope by instruction): the
  unconstrained legacy allocation constants (depth_bias_strength, direct_sulfur_bonus); the
  Ohsu kokumi EC50s, which are three orders of magnitude off their now-corrected source and
  are still pinned by a unit test as a drift guard; the README's illustrative sample-output
  table; and everything other agents surfaced after this wave started.

- 2026-08-27 (WAVE G3 — front-door UX / CLI honesty / tier vocabulary / determinism).
  Closes the "front-door smoke + determinism" findings (1)-(10) and module-hunt TIER 2
  optimizer items. No artifacts regenerated; no science moved (verified: predictions are
  byte-identical before/after).
  (1) THE HEADLINE COMMAND. Decision taken: **option (a)** — `--target` now screens the
  user's formulation as one extra candidate (`Your formulation (custom)`) alongside the
  15 grid entries via the existing `grid_override` plumbing, and says so. When `--target`
  is given WITHOUT precursors, the unusable formulation args are named in a warning that
  points at forward mode. The CLI states which candidate the report describes and where
  the user's own candidate ranked. Report §1 no longer echoes `vars(args)`: inverse design
  emits the evaluated candidate's real composition plus `candidates_screened`,
  `custom_candidate_included` and an explicit `ignored_cli_arguments` list; forward mode
  emits the formulation actually scored, with output plumbing dropped. Forward mode is now
  documented in README (new Step 1 vs Step 2) and QUICKSTART (mode-comparison table).
  (2) OPTIMIZER. `trial.params` is no longer passed in as `molar_ratios`. A shared mapper
  (`src.bayesian_optimizer.formulation_from_params` + `build_molar_ratios`) is used by BOTH
  `objective()` and the CLI, so the re-scored winner reproduces the trial (verified:
  trial target_score 0.2062 vs re-evaluation 0.21). Fixes, in one place: the flat 1:1:1
  fallback, the pH-into-L-Phenylalanine substring leak, and the unconditional
  `float('rosemary_extract')` crash. Concentrations print in **mM** (were labelled "M"),
  routed by parameter kind. Post-loop degenerate-objective detection prints an explicit
  "nothing was optimized; this is trial 0's random draw" banner when the objective has zero
  variance. `--pre-process` now says out loud that it is not applied (the optimizer searches
  pre-processing itself). `wheat_gluten` REMOVED from `--protein-type` in both CLIs: it is
  not a `ProteinType` member and killed every run in `_coerce_protein_type`; the matrix
  layer has only a placeholder for it, so it is withdrawn rather than aliased to soy
  (registry-backed `--protein-source wheat_gluten` is unaffected).
  (3) PRECURSOR SKIP. The `except ValueError: continue` now collects `{name, reason,
  unresolved_precursors, requested_precursors}` into `MaillardPipeline.last_skipped_formulations`
  and onto every `FormulationResult.skipped_formulations`; printed as a stdout block by
  run_pipeline, rendered as a warning banner in report.md, carried in report.json.
  `evaluate_single` now raises with the offending precursor named instead of a bare
  "Evaluation failed for formulation".
  (4) INTERVENTION MATCHING. `src/pre_processor.py` gained exact-token matching
  (`resolve_intervention` / `intervention_is_active`) over a `KNOWN_INTERVENTIONS` set.
  `{"yeast_fermentation": False}` and `"no_yeast_fermentation"` no longer switch it ON;
  `{"name": ...}` library-style records and `{"time_hours": ...}` params behave as before.
  (5) FORMULATION SCHEMA. `to_dict()` now emits the canonical pipeline-read keys
  (`temp`, `aw`) alongside the dataclass names, so a serialized Formulation no longer loses
  its temperature/water activity on the round trip; `"pH"` added to the input alias tuple
  and to `get()`.
  (6) QUICKSTART TRUTH. `docker_maillard.sh quickstart` used grid names that do not exist
  ("Baseline") — corrected to the real names, and `MaillardPipeline.resolve_grid_formulation`
  added (exact -> case-insensitive -> unique substring -> error listing every valid name),
  wired into run_campaign. `data/campaigns/shareable_meaty_screen.yml` CREATED (QUICKSTART
  referenced a directory that did not exist). Ingest docs now list ALL required flags —
  the audit named two, the code requires seven plus `--precursor`
  (`--source-reference` was also missing). Ingest preview mode now defaults to
  `results/ingest_previews/` instead of writing four files into `results/validation/`
  uninvited; `--confirm` keeps the old destination. The README claim that `--confirm`
  "updates the benchmark database and regenerates validation artifacts" is replaced by what
  it does (writes one YAML), stated in the README, the QUICKSTART, the `--help` text and the
  CLI's own closing line. `generate_validated_envelope_report.py` ALREADY had `--output-dir`
  (audit finding (7) was stale); the real wart — `configure_science_plot_style()` at import
  time, which broke even `--help` without dvipng — was moved into `main()`.
  (7) TIER VOCABULARY. Three vocabularies genuinely coexist and all three are now named:
  `tier` (high/medium/low/exploratory, run + per-compound + per-tag, paired 1:1 with
  `prediction_mode`), `calibration_evidence_strength` (literature_anchored ->
  conditional_literature_anchored -> class_anchored -> directional_transferred ->
  process_state_mismatch -> heuristic), and the literature-prior `confidence_tier`
  (high/medium_high/medium/medium_low/low, a property of `data/lit/` sources). The report's
  §6 glossary previously defined `bounded_calibration`/`transferred_literature`/
  `surrogate_family`/`xtb_derived` at compound level, which NO column emits — rewritten to
  the emitted vocabularies with score bands and the licensed action per tier. The
  campaign/comparison JSON key `confidence_tier` (which actually holds the run-level `tier`)
  now ships alongside a canonical `tier` key, the leaderboard column is renamed "Tier", and
  the collision is stated in the glossary. README sample-output table replaced with the REAL
  Compound Confidence table (real columns, CI embedded in the Predicted cell, real
  `exploratory` tier values) plus a note that all-`exploratory` is the honest current output
  for matrix systems. The malformed 10-field row in the 8-column table is fixed by a new
  `_md_cell()` escaper (`observable_assumption_summary` is itself pipe-joined); applied to
  the compound-confidence, evidence-ladder, calibration-summary, §1 and leaderboard tables.
  GLOSSARY.md restructured into "labels emitted in output" (~25 defined, including the 15
  that were emitted-but-undefined) and "concepts, not emitted as labels" (the 4 defined
  terms that appear nowhere in output are marked as such rather than deleted).
  (8) DETERMINISM. The hash-order-dependent summation was
  `usability_reports._fallback_sensitivity_summary`: `sum(precursor_attribution.values())`
  over a dict whose key order comes from a `set` in the recommender, so identical science
  summed in different orders (~1 ULP) under different PYTHONHASHSEED. Now summed by sorted
  key, with the driver sort tie-broken by name; the waterfall totals in reporting.py got the
  same treatment. `provenance.input_fingerprint_sha256` now hashes only the scientific
  inputs — `output_dir`/`report`/argv-shaped keys are excluded via
  `projection_utils.fingerprint_inputs`, so the same run written to two directories
  fingerprints identically (full argv is still recorded under `generator`). QUICKSTART gained
  a Reproducibility section.
  TESTS: `tests/unit/test_front_door_ux_2026_08.py` (48 tests, all passing) covering
  optimizer param mapping, degenerate-objective detection, intervention token matching,
  Formulation round-trip, precursor-skip surfacing, grid-name resolution, markdown-cell
  escaping, glossary/README truthfulness, fingerprint stability, order-independent
  summation, in-process bit-determinism, and both CLI modes' report inputs.
  Two pre-existing tests updated with dated reasons: `test_bayesian_optimizer.py`
  DummyTrial doubles now expose `.params` (the objective reads it), and the two
  `test_usability_reports.py` glossary assertions that PINNED the fictional
  `bounded_calibration`/`external_failing` vocabulary now assert the emitted one.
  KNOWN ENVIRONMENTAL: `test_cli_scripts.py::test_run_campaign_named_comparison_mode` still
  fails on the missing dvipng (unchanged, tracked above).
- 2026-08-27 (G3 addendum): campaign smoke passed end-to-end (exit 0, all artifacts).
  NEW cosmetic follow-up: the campaign leaderboard roll-up prints pH 0.00 / Temp 0.0 for
  every row — generate_campaign_report reads conditions.get("ph")/.get("temp") from
  per-run dicts that use different key spellings; per-run reports show the correct
  values, only the roll-up columns read zero. Owner: run_campaign/campaign reporting;
  candidate for Wave H or a later pass.
- 2026-08-27 (Wave G1 landed — chemistry): all 9 fix groups in. Fabricated radical
  chain removed (was 93% of the full network: 5500 -> 369 steps); hexanal/norfuraneol/
  glyoxal/CML/DMHF/correct-disulfide newly reachable; real 3-step norfuraneol MFT route
  primary (2,3-enolisation 28.0 / cyclisation 28.0 est / H2S-addition 28.60 Hofmann-
  anchored), legacy shortcut demoted + de-[S]'d; H2 fiction confined to 4 whitelisted
  lumping notes; 0 families at DEFAULT_BARRIER (test-enforced); thiamine cascade
  rebalanced; curated/engine family vocabulary unified (review's "wiring inert" claim
  corrected: defect was vocabulary split); quinone_cys prior parked; species-registry
  SMILES corrected where they contradicted their own InChI. CONSEQUENCE (the fitting
  tell resolving): sulfur absolute yields fell 5-40x — Hofmann1998 MFT 1.02x -> 7.83x,
  xylose lanes 1.2-2.3x -> 40-95x, Cerny flipped over->under. 17 scientific/integration
  failures, one root cause, deliberately NOT retuned. Carried to Wave H: Hofmann-only
  sulfur refit, snapshot refresh, regeneration, reconciliation, README 16-family table
  honesty, network-diagram colour map, MAX_MW=300 lipid cap note, curated Enolisation_1_2
  ascorbic-prior provenance [P].
- 2026-08-27 (three carried fixes landed: campaign leaderboard, network colour map,
  DOI-less identifier waivers). No science moved; no validation artifact regenerated.
  (A) CAMPAIGN LEADERBOARD ZEROES. Root cause was not only key spelling: campaign
  formulations come from data/formulation_grid.yml, which carries NO process
  conditions at all, so `conditions.get("ph", 0.0)` / `.get("temp", 0.0)` in
  `src.reporting.generate_campaign_report` could never resolve — every row printed
  pH 0.00 / Temp 0.0, and the Protein column printed "free" for the same reason while
  the run was actually pea_iso. Fix: a shared resolver
  (`src.input_normalization.resolve_condition_value` / `resolve_condition_float` over
  `CONDITION_KEY_ALIASES`) instead of scattered aliases, and the roll-up now walks the
  same fallback chain `pipeline.evaluate_all` uses — per-run override, then the
  campaign's `shared_conditions`, then the effective global `ReactionConditions`
  (passed in by run_campaign as the new `effective_conditions` argument). An
  unresolvable condition is emitted as `null` / "n/a", NOT 0.0: a reported pH of 0.00
  is a chemistry claim. Leaderboard construction extracted to
  `reporting.build_campaign_leaderboard` so it is unit-testable without the figure
  path. Verified end-to-end on data/campaigns/shareable_meaty_screen.yml (LaTeX style
  stubbed locally — dvipng still missing): all 4 rows now read pH 5.50 / 105.0 C /
  pea_iso.
  (B) REACTION-NETWORK COLOUR MAP. `scripts/generate_reaction_network.py` was keyed on
  pre-Wave-G1 family names; 7 of the 15 families `src.curated_pathways.PATHWAYS` emits
  had no colour and drew in the grey `#7f7f7f` default. Map rebuilt against the
  families the code actually emits (Enolisation_1_2 / _2_3_Amadori / _Intermediate,
  Furanone_Cyclisation, Aminoketone_Condensation, Thiol_Addition_Norfuraneol,
  Thiol_Oxidation added). Retired keys, with dated comments in the file: `Enolisation`
  (split by G1 — original colour kept by Enolisation_1_2) and `Sugar_Dehydration`
  (emitted nowhere in the repo any more — colour repointed to Furanone_Cyclisation).
  The generator now hard-fails on any uncoloured curated family and prints a coverage
  line; the edge legend lists only families with a drawn edge. Verified by running the
  generator into a throwaway ROOT: 15/15 coloured, 0 default fallthrough, 0 stale keys.
  The PNG path still needs dvipng, which is NOT installed here — nothing was installed
  and no tracked figure was regenerated.
  (C) DOI-LESS IDENTIFIER WAIVERS -> BASELINE ZERO. The four DOI-less sources now carry
  a typed `identifier` + `identifier_scheme` (+ `identifier_note` recording the
  retyping): Huang (2022) Clemson TigerPrints thesis 3936 -> `url`
  https://open.clemson.edu/all_theses/3936/; Fraser et al. (2018) -> `patent`
  US9943096B2; Cantre et al. (2007) -> `journal_locator` "Philippine Agricultural
  Scientist 90(2):143-152 (2007)" (no URL or ISSN asserted — none was verified); Liu
  (2023) PPI off-note hold-out -> `citation` (NC State thesis; no DOI or handle
  exists, and the peer-reviewed twin 10.1016/j.foodchem.2022.134998 is named in the
  note but deliberately NOT used as the bundle's anchor). `scripts/ci/citation_gate.py`
  WAIVERS is now empty, with the ratchet history kept as a comment; gate PASSes with
  0 waivers (660 DOI-bearing fields, 232 unique DOIs, was 664/236).
  SAFETY: `bench.get("source_doi")` truthiness fed three classifiers. New shared helper
  `src.benchmark_validation.matrix_source_anchor` accepts the typed pair alongside
  `source_doi`, and `_matrix_source_origin`, `_matrix_source_reference`,
  `_matrix_external_data_status` and `uncertainty_propagation._benchmark_signal_origin`
  all read it — otherwise retyping an identifier would have silently DOWNGRADED a real
  external source to "unspecified origin"/"internal_synthetic". PROOF: before/after
  classification snapshot over all 21 benchmark files (status, origin, signal_origin,
  promotable, promotion_blocker) is byte-identical; the Liu bundle stays
  `external_validation_only` (it is short-circuited by `evidence_class` anyway). The
  only delta anywhere is the human-readable `_matrix_source_reference` string for that
  bundle: "LiuThesis2023" -> the full citation. No tracked artifact embeds it.
  `scripts/ci/holdout_guard.py` PASSes (4 bundles, all three invariants).
  TESTS: new tests/unit/test_audit_remediation_carried_2026_08.py (16 tests) covering
  all three; test_benchmark_validation / test_formulation / test_uncertainty_propagation
  / test_usability_reports / test_front_door_ux_2026_08 (89) and the matrix-evidence /
  external-validation set (40) re-run green. KNOWN ENVIRONMENTAL, unchanged:
  test_cli_scripts.py::test_run_campaign_named_comparison_mode still fails on missing
  dvipng, inside generate_report's figure path, before the roll-up.
  FOUND, NOT FIXED: `src.external_validation.materialize_external_validation` rebuilds
  the hold-out bundle's provenance from the flavor-payload anchors' `doi` field, which
  the Liu anchors no longer have (they were repaired to
  `source_status: identifier_unavailable` on 2026-08-26). A regeneration would emit
  `source_doi: ""` and would not carry the typed identifier either — pre-existing
  drift, out of scope here, worth a Wave H pass that teaches the materializer the
  typed pair.

- 2026-08-27 (WAVE H — final reconciliation: sulfur refit, regeneration, test + text
  reconciliation). Five parts, all landed. This wave's headline is a NEGATIVE result,
  reported rather than tuned away: the sulfur branch cannot be recovered by any barrier
  value, and the one constant that would recover it was deliberately left alone.

  PART 1a (ROOT CAUSE, found while setting up the refit). The two knobs this wave was
  asked to refit — `thiol_addition_norfuraneol` and `furanone_cyclisation` — had EXACTLY
  ZERO derivative on every prediction. Cause: `Enolisation_2_3_Amadori`, the step Wave G1
  added to open the accepted 1-deoxyosone -> norfuraneol -> MFT route, was collecting BOTH
  the amine-nucleophile ionisation correction and the Amadori water-activity correction in
  `src/conditions.py` — by SUBSTRING MATCH ON ITS OWN NAME. That step releases the amino
  acid; the amine is a leaving group, not a nucleophile. At Hofmann's pH 5.0 / aw 0.98 the
  accidental match cost it 1e-3 * 0.62 = +6.06 kcal/mol, which kept the real MFT route
  ~1600x below the demoted one-step shortcut Wave G1 had just replaced it with. G1's
  flagship chemistry fix was inert. FIXED narrowly: new
  `_releases_rather_than_attacks_with_the_amine` exempts families that are themselves
  enolisations / eliminations / deaminations; of the ~35 families the engine emits,
  `Enolisation_2_3_Amadori` is the only one it reclassifies. Measured panel impact at the
  shipped barriers: 2 rows move (Resconi furfural -2.9%, Bolton MFT -0.03%) — the
  norfuraneol route takes over as the selected MFT path at the identical span and flux the
  shortcut had. What it restores is IDENTIFIABILITY.

  PART 1b (SULFUR REFIT, Hofmann1998 only). New reproducible generator
  `scripts/generators/refit_sulfur_barriers_hofmann.py`; record in
  results/validation/sulfur_barrier_refit_hofmann.{json,md}. Fit target: the single
  benchmark, asserted (hold-out, Internal2026, ProtocolPilot2026 all excluded by
  assertion; the xylose HVP lanes and Cerny/Bolton explicitly named as forbidden).
  Objective: mean |log10(pred/meas)| over its two rows. Three decision rules stated up
  front and applied mechanically: IDENTIFIABILITY (flat profile -> keep incumbent),
  MATERIALITY (>= 0.02 dex to move), CONSERVATIVE EDGE (adopt the LARGEST barrier within
  0.01 dex of the profile minimum, so a saturating profile cannot drag a constant to the
  bottom of its range). Result, over ranges bounded by values already in FAST_BARRIERS for
  the same mechanistic class:
    * `thiol_addition_norfuraneol` 28.60 -> 26.85 (range [23.30, 29.65]; profile min 0.6198
      dex, incumbent 0.6987; 26.85 is the conservative edge of the 23.30-26.85 indifference
      band). Its 28.60 had been inherited from `thiol_addition` and from THAT key's PRE-G1
      Hofmann window [28.10, 28.85] — a window derived when MFT was made by the shortcut,
      so it did not survive the route change.
    * `thiohemiacetal_formation` (23.3): NOT IDENTIFIABLE — exactly zero derivative over
      its whole range. Incumbent kept.
    * `furanone_cyclisation` (28.0): zero derivative BEFORE the norfuraneol move (the
      norfuraneol step at 28.60 set the route's span, so the cyclisation step at 28.0 was
      never rate-determining); after the move it becomes identifiable (span 0.069 dex) and
      the incumbent 28.0 turns out to be exactly at its own optimum — achievable gain
      0.0000 dex. Incumbent kept, now for the stronger reason.
    * `thiol_dehydration` (26.8): IMMATERIAL — 0.0096 dex available. Incumbent kept.
    * `thiol_addition` (28.60): REPORTED, NOT FITTED. After G1 it labels only the demoted
      one-step shortcut and the lumped `Thiol_Addition_H2` step; fitting it would tune the
      route this programme exists to retire.
  Objective 0.6987 -> 0.6298 dex. Hofmann MFT 7.83x -> 5.58x under, FFT 3.19x -> 3.26x.
  The script is IDEMPOTENT and the committed artifact is the convergence re-run: at the
  adopted point no knob moves again (best remaining gain 0.0099 dex), and the move itself is
  preserved in the artifact's `applied_history`.
  THE FINDING: the profiles SATURATE inside their defensible ranges. With every fitted knob
  pinned at the bottom of its range the objective is still 0.6102 dex. No barrier value
  reproduces the measured absolute yields, because the deficit is in the volatile-budget
  ALLOCATION, not in the barriers: at the Hofmann conditions the total budget is ~1050 ppb
  against a measured MFT+FFT of 542 ppb — the right order of magnitude — and furfural,
  which that experiment did not measure, takes ~78% of it.

  PART 1c (XYLOSE HVP OBSERVABILITY, re-derived not compensated). New generator
  `scripts/generators/rederive_hydrolysate_observability.py`; record in
  results/validation/hydrolysate_observability_rederivation.{json,md}. The
  `_HYDROLYSATE_SULFUR_OBSERVABILITY_PROFILES` base factors constrain the PRODUCT, so they
  were re-derived against the same two literature benchmarks they were originally fitted
  to (PMC9905368), with the barriers now moved — same method, new chemistry — rather than
  bending a barrier to compensate. Same three decision rules, plus ADMISSIBILITY: a
  surviving fraction lives in (0, 1].
    * Methional 0.0045 -> 0.05623. Interior and identified (unconstrained optimum 0.0639);
      residuals 16.2x / 12.4x under -> 1.30x under / 1.01x. Gain 1.09 dex.
    * 2-Furfurylthiol and 2-Methyl-3-furanthiol: SATURATED, NOT MOVED. Their unconstrained
      optima are 8.65 and 3.49 — 8.6x and 3.5x ABOVE the physical maximum of 1.0. A factor
      that can only suppress cannot fix an under-prediction; the layer is exhausted and the
      21-95x residual is somewhere else. Incumbents survive as unconstrained legacy
      estimates, explicitly NOT as fitted values.
    * bis(2-methyl-3-furyl) disulfide: NOT DERIVABLE — it appears in no literature benchmark
      of this lane; its only comparators are the synthetic snapshots, which are forbidden.

  PART 1d (Cerny / Bolton / acrylamide re-pins). Cerny2008 FLIPPED from over- to
  under-prediction (3.419x over -> 3.195x under); its test said "re-derive, not relax", so
  it is renamed, inverted and two-sided-pinned, with the VALUES_NEED_RE_EXAMINATION caveat
  retained and the contract tolerance untouched. Bolton1994 stays an honest failure and its
  inverted guard is re-pinned 173.3x -> 744.9x (the loose `ratio > 10` bound would have
  stopped detecting drift; it is now [300, 1500]).
  `acrylamide_spi_extrusion_130C_ACSRef3` 1.04x -> 6.42x under G4's SME desaturation: given
  a `model_change_note` in the benchmark JSON in the `unit_correction_note` style, test
  renamed and re-pinned to [4.5, 9.0]. The SME model is NOT bent back — the old 1.04x was
  read off a SATURATED surface (the 40 C offset cap sat exactly at this benchmark's
  250 kJ/kg SME), so the agreement was with a constant, not with a model. The whole 6.2x is
  exp(-(129 kJ/mol / R)(1/421.2 - 1/443.2)) on the effective-temperature drop 170 -> 148 C.

  PART 2 (REGENERATION, once). `refresh_internal_reproducibility_snapshots.py` NOTES bumped
  to v5 naming the G-wave chemistry; then benchmark_summary, benchmark_index,
  prediction_uncertainty (n=200), external_validation_report, experiment_value_ranking,
  loo_leverage, gap_heatmap + experiment_brief_cards, experiment_requests (+ 3 stale
  rank-suffixed files pruned), validation_overview and family_validation_overview md/json,
  matrix_sigma_residual_derivation, projection_constant_refit,
  family_implementation_status, hexanal_nonanal_*, objective_progress. Both figure
  generators gained `--skip-figures`, and `generate_family_validation_figures.py`'s
  import-time `configure_science_plot_style()` (which made the module unimportable without
  dvipng — even `--help`, even the md/json half) was moved into the figure path, the same
  wart Wave G3 removed elsewhere. PNGs other than gap_heatmap remain stale and the README
  now SAYS SO.
  NEW HEADLINE NUMBERS: panel 16 benchmarks / 41 matched rows / 13 in the MC panel;
  aggregate CI coverage 28/41 (68.3%, was 33/41); split = external_literature 2/11 (18.2%,
  median CI width 3.00 dex, 4 not-evaluable) vs internal_synthetic 18/18 (100%, 3.19 dex, 8
  not-evaluable); hold-out 5/8 (62.5%) = 3/3 re-scoring + 2/5 genuine extrapolation, median
  fold error 32.79x, worst 2474x; strict-ready 0/16 (was 4/16 — the panel now has none);
  out-of-CI cells 13/41; benchmarks without blocking coverage/ranking gaps 7/16 (was 11/16); worst quantitative benchmark CML/CEL at 1204x, worst reference-only Cerny at 3.195x. Read the literature
  row with its width: the intervals GREW (2.0 -> 3.0 dex) and coverage still FELL (7/11 ->
  2/11), i.e. the predictions moved further out than a ~1000x interval allows.
  PROJECTION CONSTANT — DELIBERATELY NOT APPLIED. Re-running
  `refit_projection_constants.py` on the new chemistry moves the optimum
  `reference_conversion_time_min` 1.2589e4 -> 2512 (objective 0.8777 -> 0.7543 dex). It was
  NOT applied. tau_ref is a single GLOBAL SCALE: at that optimum the Hofmann MFT residual
  collapses to 0.048 dex while Resconi furfural grows to 19x OVER. Raising the budget 5x to
  close a 5.6x ALLOCATION gap means supplying ~5250 ppb of volatiles to explain 542 ppb of
  measured product and sending the rest to the species nobody measured. The script now
  records shipped-vs-fitted, the decision and its reasoning; `src/projection.py` carries the
  same note at the constant.

  PART 3 (carried fixes). (a) acrylamide — see PART 1d. (b)(c)(d) landed in the dated
  entry above (campaign leaderboard conditions, reaction-network colour map, typed
  identifiers + citation gate ratcheted to ZERO waivers with a proven-identical
  classification snapshot). (e) NOTE-ONLY items confirmed and left unfixed: the curated
  `Enolisation_1_2` ascorbic-prior provenance mismatch was already flagged in
  `src/curated_pathways.py` (both averaged priors are ascorbic-acid priors producing
  ~6.6 kcal/mol for a step the FAST table puts at 28.0) — confirmed, untouched. `MAX_MW`
  was only a bare comment; it now carries the measurement that makes the limit visible:
  at 300 Da the cap KEEPS linoleic acid (280.5) and the alkoxy radical (295.4) but PRUNES
  the peroxyl radical (311.4), 13-HPODE (312.5), 9-HPODE (326.5) and arachidonic acid
  (304.5). So `Radical_Propagation_O2` and `Peroxy_H_Abstraction` cannot deposit their C18
  products, a hydroperoxide supplied as a PRECURSOR survives while the same molecule made
  by the network does not, and the chain propagates only as far as the alkoxy radical.
  Hexanal (100 Da) is unaffected. Not changed — it is a network-size decision for a
  dedicated lipid workstream.

  PART 4 (test reconciliation). THE THREE "PRE-EXISTING" FAILURES WERE NOT PRE-EXISTING:
  all three PASS at HEAD (5fe8515) in a clean worktree, so they were broken on this branch
  and the earlier "verified pre-existing" was relative to that agent's own starting point.
  Diagnosed properly:
    * `test_matrix_benchmark_ranking_contract` — fixed by the snapshot refresh; no longer
      failing.
    * `test_canonical_systems[ribose_cysteine_leucine]` — RE-DERIVED, not relaxed. The old
      claim was "at least one characteristic target in the top 10 by simulated concentration
      in the CANTERA lane". All five targets are still reachable and non-zero; what changed
      is their magnitude in that lane only (FFT 3.8e-7 -> 6.8e-14 mol/L). Cause: Wave G1
      removed the fabricated H2 economy, and the three steps that reach these targets all
      consume the `[HH]` REDUCING-EQUIVALENT TOKEN. The production path pool-GATES that
      token; the Cantera export takes it literally as a mass-action reactant, so with the
      fabricated H2 source gone those steps are throttled by an [H2] of ~1e-8. Also: the old
      top-10 for this system contained `MFT_radical`, `elemental_sulfur` and
      `quenched_elemental_sulfur` — it was partly pinning a ranking of fictitious chemistry.
      The test now asserts a STRICTLY STRONGER claim for all three systems (every named
      target reachable and non-zero) and keeps the abundance claim only where it is
      meaningful. Two latent defects fell out: the sanitizer never replaced commas (so
      "2,5-dimethylpyrazine" never matched) and "FFT" was never resolved, so two of the three
      systems had been passing on a single target each while their others went unchecked.
    * `test_calibration_loop_improves_mae` — a real bug, fixed. `_compute_prediction_error`
      clamped the PREDICTION side of the log error at the same 0.1 ppb floor as the
      measurement. That is right for a measurement (GC-O quantitation limit) and fatal for a
      prediction: the objective becomes exactly constant for every prediction under 0.1 ppb,
      L-BFGS-B sees a zero gradient, and `calibrate_matrix_constants` reports "did not
      improve the MAE" and reverts — silently, for a reason unrelated to the calibration.
      After G1 dropped sulfur yields that is the normal case (the observed MAE was 3.5,
      which is the value of the clamp, not of any prediction). The two floors are now
      separate named constants and the prediction floor is a pure log(0) guard at 1e-12.
  Also re-pinned with dated causes: the HVP xylose magnitudes (5.165x -> 94.6x SPI /
  46.8x wheat, ordering claim untouched), the free-AA regression (one-sided calibration
  ceilings -> two-sided pins on a quantified under-prediction, table renamed
  `BENCHMARK_EXPECTED_FOLD_ERRORS`), and the example matrix intake's expected ranks
  (2,5-dimethylpyrazine lifted above the thiols; MFT swapped ahead of FFT).
  Also re-pinned: `tests/unit/test_budget_projection.py`'s Methional observability
  assertion (0.0045 -> 0.05623), with the re-derivation named; the property that test is
  about — Methional is NOT source-sensitive while the thiols are — is unchanged.
  NEW: tests/unit/test_wave_h_2026_08.py (9 tests) pins the ionisation exemption (both that
  it fires on the enolisation and that it does NOT fire on the genuine amine-nucleophile
  steps), the calibration-objective gradient, the two fits against their records, the
  (0, 1] bound on every observability factor, and that both fit scripts REFUSE synthetic and
  hold-out targets, and that the README's "5 of 16" claim still matches the derived
  artifact (so the family table cannot drift the way "16 families of reaction chemistry"
  did). Both fit scripts are IDEMPOTENT and carry an `applied_history` block, so the
  committed artifact is the convergence re-run while the move it made stays on the record.
  FINAL SUITE (certifying run on the settled tree): 1178 passed, 5 skipped, 2 xfailed,
  1 failed — the failure being
  `test_cli_scripts.py::test_run_campaign_named_comparison_mode` on the missing dvipng,
  unchanged and environmental.

  PART 5 (text coherence + the 16-family honesty). New generator
  `scripts/generators/generate_family_implementation_status.py` derives the answer from the
  engine by enumeration over a battery of entry points, asserts that every emitted reaction
  family is classified (a new template with no mapping fails the script), and writes
  results/validation/family_implementation_status.{json,md}. RESULT: **5 of 16** lanes carry
  generative reaction templates — 01 amino acid + sugar core (20), 02 lipid oxidation (5),
  03 thiamine (2), 11 lipid x Maillard crosstalk (3), 12 protein damage markers (4). The
  other eleven are shared limbs (07, 08, 09), matrix/process layers that are not reaction
  chemistry (06, 10, 16), the curated layer only (14), or literature priors with no
  generative implementation (04, 05, 13, 15). README's family table gained an
  Implementation column saying exactly that, the "16 families of reaction chemistry" prose
  claims and the mermaid node were corrected, and two entry-point facts were recorded that
  the table cannot show: the lipid radical chain enumerates to ZERO steps from an unoxidised
  fatty acid plus O2 (it needs a hydroperoxide seed, so in production the anchor seeds it),
  and the thiamine cascade matches thiamine by EXACT canonical SMILES plus a >= 100 C gate.
  README calibration badge/headline/split table, the surfaces table, the trust table (the
  "High" tier is GONE — free-precursor sulfur chemistry used to be it) and a figure-staleness
  note updated; VALIDATION_CONTRACT §2 (strict-ready 0/16, with the 6 -> 4 -> 0 sequence and
  the note that no tolerance was widened) and §3F (current split, with the read-it-with-the-
  width warning); CURATING_LIPID_OXIDATION_ANCHOR gained a dated note on the radical-chain
  rebuild and the hydroperoxide-seed requirement.

  [P] CARRIED FORWARD, not addressed here: `direct_sulfur_bonus` in `src/recommend.py` is
  now DEAD — it tests `terminal_family == "thiol_addition"` and Wave G1 renamed every
  terminal sulfur family, so an unconstrained legacy heuristic switched itself off in the
  same commit that changed the chemistry. Measured size: at most 1.68x, and 1.007x at the
  Hofmann conditions, so it is NOT the cause of the sulfur collapse. Left dead ON PURPOSE
  and documented in place: it is the exact knob that would absorb the MFT residuals this
  panel exists to expose, so re-pointing it belongs to the allocation workstream, with a
  refit. That workstream (`depth_bias_strength`, `direct_sulfur_bonus`, and the
  budget/allocation split) is now the highest-value open item: it is what the sulfur refit
  identified as the location of the residual. Also carried: the
  `materialize_external_validation` typed-identifier drift flagged in the entry above; the
  Ohsu kokumi EC50s; the CML/CEL and furosine unit-collision gaps (1204x / 201x), which are
  real predictor gaps rather than reporting artifacts.
- 2026-08-27 (CLOSE-OUT): Wave H landed (refit honest: barriers cannot fix the sulfur
  deficit — allocation-layer residual documented, papering-over knob declined; final
  numbers 2/11 lit rows in CI / 0/16 strict-ready / 5/16 generative lanes; suite 1178
  passed / 1 dvipng; gates 0 waivers). AUDIT.md written (public narrative) + README
  audit banner. Work committed as 7 thematic commits on audit-remediation, pushed,
  PR #12 opened against s27 for self-review (scratch/ + .claude/ excluded and
  gitignored). Cold-start red team launched against the committed tree: a zero-context
  scientific due-diligence reviewer and a zero-context forensic code auditor — the
  final verification that the remediation holds under the same class of scrutiny that
  started it.
- 2026-08-27 (RED TEAM, forensic report — the step worked): verdict "not audit-washing,
  but incomplete in the direction that flatters the model; honest count 0/11".
  FABRICATION-GRADE: (F1) the 2/11 headline's only hits are the two methional rows the
  hydrolysate observability constant was FITTED to during Wave H (a circularity my own
  wave instruction introduced — re-derive-and-then-score); (F2) both PMC9905368
  benchmarks describe an experiment the cited paper never ran (fructose/glucose pH 7.5
  relative-peak-area paper vs our xylose/cys pH 6.0 absolute-ppb files) — the citation
  sweep's identity-not-content blind spot, under 6/15 literature rows; (F3) Trikusuma
  8-sig-fig undisclosed back-fit (soy disclosure also missing); (F4) Pratap-Singh
  hexanal ~4.35x vs paper Table 1 (reported, unconfirmed — MDPI 403). MISLEADING: stale
  README 0.26 dex (artifact says 0.75); bi_2020 hold-out 1260 ppb is REPO-COMPUTED from
  OAV x our own threshold; family_sensitivity false-zero for schiff/enolisation/cys
  (offset key routing bug); MC clamp binds at sigma 2.86 (+-100x actual vs +-110x
  advertised, 10.9% draws pinned); lipid saturation docstrings claim an enabled cap
  that ships null; legacy MFT shortcut still fires for pentoses; no_verifiable_source
  population is 93 not 43, 48 numbers still act; citation-gate line overstates what the
  blocking gate can test; guardrail-test subprocess validates baseline not candidate.
  HELD UP: hold-out partition (leak hunt failed), bit-identical reproducibility,
  quarantine forensics "exemplary", ordering claims sound, unit chains coherent.
  -> WAVE I fix list queued; awaiting the scientific due-diligence red team report.
- 2026-08-27 (RED TEAM, scientific due-diligence report): FUND-WITH-CONDITIONS (~65%).
  Converged independently on F1-F4 (PMC9905368 fabrication C1, methional self-fit C2,
  back-solved anchors C3, synthesized hold-out points C4, saturation-cap headline
  preservation C5). New: (H1) flagship mechanism DOI 10.1021/jf60200a038 resolves to a
  gossypol/rat paper — correct is jf60199a045; gate blind to src/docs/README; (H2)
  pentose>>hexose is 7.05x not the published 15.8x, and carried by the unconstrained
  hexose legacy fit (zeroing the 2.8 kcal gap -> 1.69x), not structure; (H3)
  Lipid_Schiff_Base SMARTS = 83% of the core network (exclusion comment not
  implemented); (H4) MFT is gated on pool-H2 produced only by pyrazine chemistry;
  (H5) broken Arrhenius pairs (amadori t1/2 0.13 ms; thiol_addition faster than
  diffusion; Schiff/Amadori ordering INVERTS between FAST_BARRIERS and the YAML;
  lipid_homolysis kinetically dead); (H6) off-note flagship non-functional in the user
  path (no fatty-acid precursor; --lipids only accepts the aldehydes themselves;
  quickstart reports Off-Flavour Risk 0.00); (H7) MFT:FFT ordering inverted vs
  measurement on the hydrolysate rows; (M1) origin/main still shows brightgreen 39/48
  + no AUDIT.md; (M2) README showcase table unreproducible + stale May "TRUSTED"
  example artifact; (M4) VoI is a heuristic, not value-of-information; (M5) rank-1
  experiment card prescribes the wrong protocol (substring-matched template); (M6)
  --report hard-crashes without dvipng on the documented conda path; (M7) count
  mismatches (76-papers claim vs 42 entries; no_verifiable_source 43 vs 66-93);
  (M8) "Safety Score: 1.04" unlabelled (2x risk, higher=worse). Conditions + maintainer
  questions logged in the red-team transcript. -> WAVE I launched.
- 2026-08-27 (due-diligence ADDENDUM): reviewer re-ran the full suite on Docker (has
  dvipng): 1202 passed / 2 failed — both stale assertions pinned to the just-quarantined
  PMC9905368 files; withdrew any implication the test claim was inflated. Observed the
  Wave-I remediation live and judged the self-incriminating quarantine write-up "the
  single strongest piece of evidence I gathered" — verdict nudged toward the optimistic
  end of fund-with-conditions. NEW CONDITION 8: land Wave I as ONE atomic green commit
  (quarantine + regeneration + republished 0/9 headline + badge together). NEW
  MAINTAINER QUESTIONS: (Q8) what blocking check replaces identity-only citation
  verification; (Q9) what is the prior that the remaining 9 literature rows are clean
  and what would raise it (candidate answer for close-out: add a per-benchmark
  content_verification status — identity_verified vs content_verified-against-fulltext —
  and a campaign to full-text-verify the 9 survivors); (Q10) what prevents test re-pins
  from re-encoding sourceless numbers.
- 2026-08-27 (RED TEAM, test-suite integrity report): ~15-20% of suite is real
  assurance, ~2% scientific; after 66c0ace no test asserts out-of-sample quantitative
  agreement (all converted to honest pinned-miss bands — verified no contract widened).
  CRITICAL: (S1) "1178 passed" NOT reproducible from a clean checkout — 6 failures on
  the exported certified commit because .gitignore data/* excludes required test data
  (data/qm/*.json, data/qm/irc_validation_cases/, data/ingest_templates/*.csv,
  data/Gemini_Deep_Research/); the certified number was local-only; (S2) CI HAS NEVER
  PASSED — unit tier fails on the missing data, tests-full is gated behind it so the
  scientific/integration tier has NEVER run in CI; PR #12 is red; (S3)
  test_benchmark_correlation computes no correlation; MAILLARD_STRICT_BENCHMARKS set
  nowhere in CI/Make; (S4) test_trikusuma..._is_quantitatively_supported asserts
  <=1.05x against 3 same-commit back-fitted constants (fabrication-adjacent name);
  (S5) Pouvreau-2021 pH test: source uncited anywhere in repo, pins a tuned knob
  (exp(2x0.235)=1.600 dead-centre of the band), pyrazine assertion tautological;
  (S6) data/qm barrier windows entirely uncited + outside citation-gate scope, echoed
  by conftest fixture "literature_barriers_dict"; (S7) 15 pass-body tests
  (test_cantera_sim), 5 no-op balance parametrizations, 2 self-excusing skips, 2 dead
  imperative xfails; 19.4% of tests all-trivial; 36% of tests/scientific is
  report-shape only; old assert_balanced drops unparseable species (10 vacuous
  mass-conservation tests); (S8) tests/benchmarks QM authority lane deleted Apr 21,
  orphaned scaffolding + loader-echo tests failing on unshipped data; (S9) NOTHING
  test-guards hold-out accuracy or the headline numbers (badge can drift silently).
  CREDIT: re-pins verified honest; invariant core "genuinely excellent";
  pentose>>hexose test "the model of how to do this". -> WAVE J queued (after Wave I
  lands): ship the missing test data (targeted .gitignore un-ignores), make CI
  actually green end-to-end, rename/fix the deceptive tests, delete/revive pass-body
  tests, fix old assert_balanced, add headline+hold-out guard tests, cite-or-flag
  data/qm windows, extend gate scope to tests/ and data/qm/.

- 2026-08-27 (WAVE I — remediation of the cold-start red-team findings). Structure:
  fixes 1-5 CRITICAL, 6-12 HIGH, 13-18 MEDIUM/DOC, then regeneration + reconciliation.
  Two parallel sub-agents (chemistry/DOI stream; tooling/bug stream) with strict file
  ownership; everything else in the main session. FULL PER-FIX RECORD IS APPENDED BELOW
  ONCE THE REGENERATION LANDS -- this entry is the structural record.

  FIX 1 (QUARANTINE + REVERT). `spi_hvp_xylose_120C_PMC9905368` and
  `wheat_gluten_hvp_xylose_120C_PMC9905368` git mv'd to data/benchmarks/quarantined/ with
  a full evidence block. THE FINDING: their DOI is LIVE and its metadata MATCHES -- which
  is why the 225-anchor sweep passed them -- but the paper (10.1007/s10068-022-01194-w)
  reacts hydrolysates with GLUCOSE and FRUCTOSE at pH 7.5 for 90 min and reports only
  RELATIVE PEAK AREAS, and never mentions FFT or MFT. A relative-peak-area paper cannot be
  the source of an absolute ppb value by any route. Six numbers with no possible
  provenance. Methional base_factor REVERTED 0.05623 -> 0.0045 in
  src/recommend.py::_HYDROLYSATE_SULFUR_OBSERVABILITY_PROFILES; the whole table is now
  labelled as having NO surviving literature constraint of any kind.
  results/validation/hydrolysate_observability_rederivation.{json,md} marked RETRACTED
  (banner + JSON block) and its generator refuses to run (exit 2, named reason).
  projection_constant_refit gained a `contamination_note` naming the removed rows (6 of 15)
  and what survives. sulfur_barrier_refit_hofmann gained a `contamination_review`:
  VERDICT UNCONTAMINATED, RESULT STANDS -- it fits cys_ribose_140C_Hofmann1998 only
  (10.1021/jf9705983, verified) and the quarantined files were already on its forbidden
  list, contributing zero rows. Standing caveat recorded: Hofmann1998's own 1.45/0.09
  contract is the 3rd-tightest in the collection, and it is now the ONLY literature anchor
  on the entire sulfur branch.

  FIX 2 (BACK-SOLVED CONSTANTS). New evidence strength `fitted_to_benchmark` in
  src/matrix_calibration_registry.py, applied to 11 records across three lanes, each with
  the arithmetic written out so the claim is checkable rather than asserted:
  pea reference yields x 1268-1271 = that benchmark's own measured 260/638/80 ppb;
  soy numerators x 838.5-839.2 = 380/2492/120 ppb; Trikusuma's three 17-sig-fig constants
  score max_ratio 1.000 / MALE 0.000 on their own benchmark. 17-sig-fig fit residue rounded
  to 6 sig figs with `previous_value` retained (relative change < 1e-6). New
  `BenchmarkSummary.evidence_role` (predictive | fit_recovery | internal_synthetic); the
  summary markdown gains an Evidence Role column, marks those rows
  "pass (fit-recovery, not predictive)", and the "N/16 without blocking gaps" line is
  replaced by a 4-line breakdown whose headline is the PREDICTIVE count.

  FIX 3 (HOLD-OUT HONESTY). Every hold-out row now carries `value_provenance` +
  `value_provenance_note`: `band_geometric_midpoint` (the Liu rows -- sqrt(min*max) of a
  10-12x band, a construction of ours) or `derived_from_oav_and_repo_threshold` (the Bi
  rows -- the source's OAV x THIS REPO'S OWN unsourced 4.5 ppb hexanal ODT; 280 x 4.5 =
  1260.0 and 72 x 4.5 = 324.0, to the digit). Values rounded to 4 sig figs;
  band_min/max_ppb carried; `uncertainty_pct` kept for schema compatibility but labelled as
  the UPWARD half-width of a geometric band, not a symmetric analytical uncertainty. The
  fabricated `measurement_date` ("<publication year>-01-01") is gone -- `not_applicable`
  plus a note saying what used to be there; `_publication_year_proxy` retired in place with
  the reason. The external report now leads with coverage at the PRE-WIDENING sigma (2.0)
  computed by re-running the identical hold-out through a new
  `default_priors(matrix_sigma_override=...)`, states the shipped-sigma figure as a
  secondary prior-dependent line, and renders a measurement-vs-derived split table.

  FIX 4 (FLAGSHIP DOI) + FIX 5 (NEW GATES): see the sub-agent record and
  scripts/ci/fit_target_gate.py. The gate is negative-tested (both defect classes injected
  and caught) and wired into ci.yml's dependency-free tier as BLOCKING.

  NEW MACHINERY: src/fit_target_index.py reads the fit records under results/validation/
  and classifies each by LEVERAGE (free parameters per fitted row). `per_row_recovery`
  (>= 0.5) excludes a benchmark's rows from the literature-coverage numerator AND
  denominator; `global_low_leverage` does NOT, because excluding a one-global-constant fit
  across two dozen rows would delete genuine FAILURES rather than expose them. The rule is
  deliberately symmetric about which direction it flatters. Undeclared leverage is a gate
  failure, not a default. Consequence, stated rather than smoothed: with the projection
  scale fitted across the whole literature panel, NO literature benchmark is fully
  out-of-sample, and the only genuinely out-of-sample evidence is the hold-out set.

  FOUND WHILE FIXING (neither red team saw these):
  (a) 10 of 14 Monte-Carlo BARRIER PRIORS WERE EXACT NO-OPS -- `schiff_condensation`,
      `1,2-enolisation` (which could never have matched anything: its target contains a `-`
      that key normalisation has already replaced), `2,3-enolisation`,
      `cysteine_thermolysis`, `thiol_addition`, `thiol_addition_hexose`, `retro_aldol`,
      `dehydration`, `beta_elimination`, `lipid_thiazole`. Same root cause as the
      family_sensitivity false zeros: `get_barrier()` resolves BARRIER_OFFSETS against the
      normalised RAW label before canonicalising. ~70% of the MC barrier channel was inert
      and every published CI was narrower than the priors it claimed to sample. FIXED via
      `_routed_offset_keys` (reuses family_sensitivity's resolver so the two cannot drift).
      MEASURED: median CI width 0.993 -> 1.155 dex, mean 1.32 -> 1.72 dex; aggregate
      coverage 27/35 -> 29/35 at n=60. NOTE THE DIRECTION: this made coverage too LOW, so
      the fix RAISES it -- reported as an interval-width change, never as accuracy. No
      prediction moved. `dehydration` remains inert ON PURPOSE and is named in a test: the
      canonical family exists in FAST_BARRIERS but the engine emits `Sugar_Dehydration` /
      `Thiol_Dehydration`, which canonicalise to themselves. Which one the MC should
      perturb is an owner decision that changes every interval.
  (b) The same routing defect in src/refinement_campaign.py's cheap-screening loop (fixed;
      it was screening five families it could not perturb and reporting the resulting zeros
      as insensitivity). [P] the ROOT CAUSE in src/barrier_constants.py (offset lookup
      before canonicalisation) is untouched -- both consumers are fixed at the producing
      side, but any future caller using a short key is still silently broken.
  (c) REGENERATING THE HOLD-OUT BUNDLES DELETED THEIR OWN PROVENANCE. The materializer
      emitted a fixed key set, so a regeneration dropped the Liu thesis's typed
      identifier/identifier_scheme/identifier_note and wrote `source_doi: ""` in its place
      -- which would have reclassified a real external hold-out as internal_synthetic --
      and erased the hand-written `citation_audit_note` on the Bi bundles. This was the
      "FOUND, NOT FIXED" item at the end of Wave G3; it is fixed now (provenance
      carry-through + `provenance_overrides` in the bundle SPEC so a regeneration cannot
      delete them).
  (d) A MIS-APPLIED DOI REPAIR, exposed by (c). The 2026-08-26 sweep wrote an EGCG /
      deoxyosone-trapping repair (basis: "EGCG lowers the activation energy for
      xylose/alanine -> ARP conversion from 77.8 to 62.8 kJ/mol") onto the two BI 2020
      RAW/ROASTED PEA HEXANAL anchors. The contradiction is inside the same record: its own
      citation_caveat says 10.1021/acs.jafc.9b07711 "should stay attached only to the pea
      hexanal baseline". Reversed with a `doi_repair_reversal` record; the superseded
      records are retained.
  (e) A CONFABULATED DOI IN docs/, surfaced the moment fix 4 widened the citation gate's
      scan to prose: `Hao et al. (2025)`, an impossible zero-padded Elsevier article
      number, scored 4/8 as a CALIBRATION source. FIXED AT SOURCE rather than waived: the
      entry is WITHDRAWN in the document (score struck, "directly calibrates the optimizer"
      claim retracted -- verified no constant is fitted to it), and the gate's
      TEXT_SURFACE_WAIVERS is back to ZERO.

  --- WAVE I CLOSE-OUT (2026-08-27) ---

  SUB-AGENT STREAMS (full per-fix records in their reports; summarised here).
  CHEMISTRY/DOI STREAM (fixes 4, 8, 9, 10, 12, 18):
    * FIX 4: 10.1021/jf60200a038 -> 10.1021/jf60199a045 in all 5 owned sites with a canonical
      doi_repair record. citation_gate gained a TEXT surface (src/**/*.py, docs/**/*.md,
      README.md, AUDIT.md) running the syntax + confabulation checks. It found ONE real
      defect immediately -- a confabulated Elsevier DOI in docs/ -- which the main session
      then fixed at source, so both waiver lists are EMPTY. Note the gate's own blind spot,
      stated by the agent: the bad flagship DOI was syntactically valid and carried no
      confabulation signature, so the widened gate could not have caught it either.
    * FIX 8 (the biggest result of the wave): the norfuraneol + H2S -> MFT step drew its 2[H]
      from a pool whose only source was pyrazine chemistry. Sourced from the thiol redox
      couple instead (2 cysteine -> cystine + 2[H], exactly atom-balanced -- dropping the
      token would have manufactured 2 H from nothing, the move Wave G1 removed). With
      aminoketone condensation disabled, MFT 0.00 -> 1038.80 ppb. On Hofmann1998, NO BARRIER
      CHANGED: MFT 61.25 -> 235.32 ppb vs 342 (5.58x under -> 1.45x under); FFT 61.44 ->
      219.96 vs 200 (3.26x under -> 1.10x OVER).
    * FIX 9: Lipid_Schiff_Base finally implements its comment's exclusions. Core no-lipid
      network share 75.2% (228/303) -> 27.2% (28/103). RESIDUAL REPORTED, NOT HIDDEN: all 28
      survivors come from one substrate, 5-aminopentanal, an omega-amino aldehyde that the
      alpha-carbon SMARTS cannot exclude; needs a whole-molecule condition. [P]
    * FIX 10: loud audit flags, ZERO runtime values changed (test-enforced). Half-lives
      computed at 150 C: amadori t1/2 = 0.133 ms; beta_elimination_dha 0.72 ms AND its
      A_value is `.nan` (silently substituted 6.25e12 -- the real defect is the prefactor,
      not the barrier); lipid_homolysis t1/2 = 1086 YEARS. The red team's "thiol_addition is
      faster than diffusion" did NOT survive arithmetic (1.6e9 L/mol/s is 1-2 orders BELOW
      the aqueous ceiling) -- the weaker true claim is written instead, plus a real defect
      the red team missed (its A_unit is bimolecular, its NaN substitute unimolecular).
      Schiff/Amadori ordering inversion flagged in BOTH tables with authority stated
      (FAST = screening lane, YAML = Cantera lane) and marked an open owner item: the two
      disagree by ~6.6e8 in the ratio.
    * FIX 12: legacy shortcut restricted to hexose (not deleted -- deleting it makes MFT
      unreachable from glucose+cysteine). Pentose>>hexose ratio unchanged: the lump never
      won the flux comparison.
    * FIX 18: template now consumes the sugar, exactly balanced. HONEST NEGATIVE RESULT:
      predicted acrylamide is unchanged to 12 s.f. (recommend.py scores from reactants and
      the sugar was already one), and SUGAR-IDENTITY dependence is still absent in BOTH
      layers -- adding one would be an unanchored new calibration. [P]
  TOOLING STREAM (fixes 6, 16, 17):
    * FIX 6: `schiff`/`enol`/`cys` (+ `retro_aldol`, `thiol_addition`, found by the agent)
      were exact no-ops. Controlled A/B on one tree: 1,2-enolisation 0.00 -> 8.86 (now the
      TOP family in the sweep), schiff 0.00 -> 0.01 (max |MAE shift| 1.58 ppb),
      cysteine_thermolysis 0.00 -> 0.00 -- but now a REAL zero, measured, not a routing
      artifact.
    * FIX 16: `configure_science_plot_style` falls back to mathtext with a one-time warning;
      MAILLARD_STRICT_LATEX=1 restores the raising contract. The dvipng test PASSES.
    * FIX 17: the guardrail subprocess re-imported from disk and so validated the BASELINE.
      Candidate constants are now passed through a shim before pytest imports anything, and
      a second latent bug was fixed (`minimize` leaves globals at its last evaluated vector,
      not res.x). Restore is now try/finally on ALL paths. The regression test was verified
      discriminating against the OLD code path.

  MAIN SESSION (fixes 1, 2, 3, 5, 7, 11, 13, 14, 15 + the four found-while-fixing items).
    * FIX 7: clamp derived from sigma at 3 sigma. Truncation 10.74% -> 0.27%; the calibrated
      tier is untouched (its clamp is unchanged at its sigmas).
    * FIX 11: [P] DELIBERATELY NOT REGISTERING a fatty-acid precursor, with the measurement
      that justifies it: the radical chain enumerates to ZERO steps from an unoxidised fatty
      acid plus O2, and MAX_MW=300 prunes 13-HPODE (312.4) and the peroxyl radical (311.4),
      so a registered precursor would feed a chain that cannot propagate and would report a
      CONFIDENT ZERO. README + QUICKSTART now say plainly what the off-note lane can and
      cannot do.
    * FIX 14: generate_report_visual_examples.py REWRITTEN to drive the real pipeline. It
      had built FormulationResults from hard-coded numbers and never touched the model --
      announcing "HIGH_CONFIDENCE (82/100)" and "TRUSTED" against a README saying there is
      no high tier, and citing a sulfur-free matrix anchor as the calibration source for the
      two flagship thiols. The real run says EXPLORATORY / tier `low` (47.0). Also fixed:
      the two README showcase PNGs were UNTRACKED (gitignored by `*.png`) and therefore
      broken on a fresh clone -- now force-added.
    * FIX 15: template selection now reads the benchmark's SYSTEM, not just the compound
      name. New `free_precursor_sulfur_yield` DoE template. The rank-1 card no longer
      prescribes "Matrix (SPI, PPI)" for an aqueous free-precursor system, and its
      dynamic-range line now spans BOTH the published 13 ppb and the model's 0.018 ppb,
      naming the 724x disagreement as the reason the experiment is ranked.

  FINAL NUMBERS (Wave I regeneration; every one of these is WORSE than what this repo
  advertised yesterday, and each fall is a fabricated support being removed):
    * Panel 14 benchmarks (was 16). MC panel 11 benchmarks / 35 matched rows.
    * Coverage split -- external literature 1/3 (33.3%, median CI width 0.86 dex, 4 not
      evaluable) | FITTED ROWS 2/2 (excluded from both numerator and denominator; BOTH
      would previously have been counted as literature hits) | internal synthetic 18/18.
      The aggregate 29/35 is retained but labelled do-not-quote.
    * Benchmarks without blocking gaps: 0/6 PREDICTIVE. (3/4 fit-recovery, 4/4 synthetic --
      i.e. every "pass" in the panel sits in a non-evidence bucket.) Was reported as 7/16.
    * Strict-ready 0/14.
    * Hold-out: 0/5 genuine extrapolations at the PRE-WIDENING prior; 2/5 under the shipped
      one with the SAME predictions. Only 4 of the 8 points are measurements at all.
      Median fold error 32.79x, worst 2474x.
    * Pentose>>hexose MFT 8.98x (was published as 15.8x) -- and only ~1.29x of it is
      structural: zeroing the 2.8 kcal/mol hexose-vs-pentose barrier gap collapses it, so
      ~7x is carried by an unconstrained legacy fit.
    * no_verifiable_source: 84 (was published as 43), of which 62 still carry numeric
      payloads consumed at runtime.
    * slr_benchmark_evaluation.md: 41 paper entries / 37 scored / 27 distinct DOIs. The
      "76 peer-reviewed papers" claim is WITHDRAWN, not adjusted -- it has no support.

  SUITE: 1245 passed, 2 skipped, 2 xfailed, 0 FAILED (tests/unit + scientific + integration
  + scripts). The dvipng failure that had been carried as "known environmental" for the
  whole audit is GONE -- fix 16 made the documented conda path work, and every figure
  generator now runs on it.
  GATES: citation_gate PASS (74 files, 846 DOI fields, 288 unique DOIs, ZERO waivers on both
  lists), holdout_guard PASS, fit_target_gate PASS (new, blocking, negative-tested).

  [P] OPEN, CARRIED FORWARD:
    1. RE-RUN THE SULFUR REFIT. Wave H's "no barrier value can fix this; the residual is an
       allocation deficit" was measured against the H2-throttled network. The saturation was
       an artifact. `thiol_addition_norfuraneol` ships at 26.85, fitted against a network
       that no longer exists. The re-run record is written and shows a further move
       available (furanone_cyclisation 28.0 -> 26.65, thiol_dehydration 26.8 -> 27.7,
       objective 0.1018 -> 0.0558 dex) and is explicitly NOT APPLIED: stacking a refit on a
       chemistry change in one pass makes it impossible to attribute the agreement.
    2. Reconcile the Schiff/Amadori ordering inversion between FAST_BARRIERS and the YAML.
    3. Lipid_Schiff_Base omega-amino-aldehyde oligomerisation (the residual 27%).
    4. Real fitted prefactors for `thiol_addition` and `beta_elimination_dha` (both `.nan`).
    5. `lipid_homolysis` 42 kcal vs the separate empirical lipid lane.
    6. The `dehydration` MC prior is inert BY DESIGN and named in a test: the engine emits
       Sugar_Dehydration / Thiol_Dehydration, which canonicalise to themselves. Which one
       the MC should perturb changes every interval and is an owner call.
    7. The ROOT CAUSE of the routing family of bugs is still live in
       src/barrier_constants.py (offset lookup before canonicalisation). Both known
       consumers are fixed at the producing side; a future caller using a short key is
       still silently broken.
    8. Reversible Schiff-base chemistry in the Cantera export (3-methylbutanal is driven to
       the numerical floor by six irreversible sinks; the concentration claim is suspended
       for that one target, named and reasoned, while the stronger structural claim -- some
       step must produce it -- is now asserted for every target with no exemptions).
    9. The lipid aldehyde-conversion cap remains DISABLED. Its rationale ends "so the
       headline stays 37/48"; that clause is now quoted in the README rather than
       paraphrased away. Re-deciding it on the evidence alone is an owner item.
   10. Sugar-identity dependence in acrylamide (absent in both layers).
   11. `direct_sulfur_bonus` / `depth_bias_strength` -- unchanged, still dead/unconstrained.

## Wave K — benchmark CONTENT verification (2026-08-27, read-only wave)

Full report: results/validation/content_verification_wave_k.md (verbatim copy of the
wave's output; repo data untouched by the wave itself — fixes queued for Wave M).
Method: every live-panel + hold-out benchmark checked against the actual paper text,
not just DOI identity. Europe PMC fullTextXML cracked all four MDPI sources that 403
direct access. 4 fatal + 2 serious findings on ~24 checkable claims:

1. FATAL — Pratap-Singh (Molecules 2021, 26, 4104, Table 1, PMC8271896): hexanal is
   1138.00±297.30 ppb (pea) and 1621.71±159.69 ppb (soy); repo uses 260/380 —
   4.38x/4.27x LOW. Independently re-verified by orchestrator via Europe PMC.
   Not a units error (same files' 2-pentylfuran matches the paper exactly).
   CONSEQUENCE: matrix reference lane anchored ~4.3x low; projection_constant_refit.md
   fitted to the wrong numbers (residual ~0.001 dex = fit recovery of wrong data);
   the S27 "matrix path over-predicts 36x" framing must be re-derived.
2. FATAL — same files: hexanol 80/120 ppb is FABRICATED. Paper Table 1: n.d. for both;
   text: "pea proteins contained no alcohol compounds"; soy total alcohols 40±9 ppb
   (1-octen-3-ol only). Both files rank hexanol expected_rank 3.
3. FATAL — Li 2026 hold-out: two of four points are WRONG TABLE ROWS (off-by-one):
   nonanal 29.42 is the Decanal row (true nonanal 74.37); 2-pentylfuran 221.5 is the
   Maltol row (true 2-PF 5625.80, i.e. 25x low). Labeled "reported_point_value".
   Hexanal + 1-hexanol in the same files are correct.
4. FATAL — Ma 2024 acrylamide 62.62 ppb not in the paper (Figure 2D only, ~150 ug/kg
   at 130C). Every plausible derivation from companion Fu 2023 ruled out. File id
   says ACSRef3 against an MDPI DOI — value predates the DOI, never re-checked.
5. SERIOUS — Hernandez furfural 1040 ppb = mean of 2 of 3 products, silently dropping
   Impossible Burger (64.71); honest 3-product mean 715.2. Cherry-pick inside +/-5% tol.
6. SERIOUS — Hofmann & Schieberle 1998 abstract reports mol %, repo uses 342/200 ppb
   with no documented conversion, on the panel's tightest contract.
CLEAN: Bi 2020 hold-out DOI correctly repaired (real pea paper; "hexanal OAV=280"
verbatim in abstract); its 1260 ug/kg point remains 50/50 OAV-derived vs measured.
UNVERIFIABLE (need library): Trikusuma 2020 (last unchecked pillar of matrix lane),
Bi 2020 tables, Liu 2023, Bolton 1994, Ramirez-Jimenez 2000.
BASE RATE: ~25% content-error rate on checkable claims -> the ~200 unchecked data/lit
anchors must be presumed similarly contaminated until content-verified.

QUEUED FOR WAVE M (after J1/J2 land, before atomic commit): correct Pratap-Singh
hexanal to 1138/1621 (+ ranking contracts to n.d. hexanol), fix Li 2026 wrong rows
(74.37 / 5625.80), re-provenance Ma 2024 acrylamide as figure_readoff ~150 with
widened tolerance or quarantine, recompute Hernandez furfural honest mean, flag
Hofmann mol%->ppb conversion gap, then ONE regeneration + doc/guard-test sync.
[P] the matrix-lane/projection re-derivation on corrected anchors is a science
decision (the old fit was recovering wrong data).

  --- WAVE J2 (2026-08-27): THE TEST SUITE ITSELF ---
  Scope: tests/ exclusively (+ two files outside it, both named below). Premise: a green
  suite that tests nothing launders confidence, so the audit's own instrument was audited.

  0. PREDECESSOR'S WORK (unstaged, unlogged) — REVIEWED, RUN, KEPT IN FULL.
     tests/unit/test_smirks_balance.py + tests/unit/test_smirks_engine.py. Verified
     empirically: 61 passed. It is correct and the reasoning in its comments checks out.
     What it does: (a) `assert_balanced` no longer `continue`s past species RDKit cannot
     parse and no longer keeps a hydrogen-suppressed molecule when AddHs raises — both
     dropped atoms from one side of the ledger, so an unbalanced reaction with an
     unparseable participant PASSED. Failures are now collected and raised naming species,
     side and family, BEFORE the element comparison (the two defects need different fixes).
     (b) New `assert_all_balanced` adds an emptiness check: ten `test_mass_conservation`
     tests were `for step in self.steps: assert_balanced(step)`, vacuously true if
     enumeration returned nothing. (c) Three negative tests pin the hardened helper so the
     laundering path cannot return. (d) The five-case `test_smirks_rule_balance`
     parametrization — the red team's "no-op balance parametrizations" — built `rxn` and
     then executed a bare `pass`. Rewritten as
     `test_smirks_rule_template_conserves_its_mapped_atoms` (parse + RDKit Validate + atom
     maps conserved, not duplicated, no unmapped heavy atoms invented), renamed because it
     is NOT full stoichiometric balance. (e) `if step:` in the Strecker test and missing
     emptiness guards in the Amadori/beta-elimination tests tightened.
     TASK 4 RESULT — NO CHEMISTRY DEFECT SURFACED: 329 steps across the 12 enumeration
     scenarios, ZERO unparseable species, ZERO AddHs failures, ZERO unbalanced steps. The
     hardening removes a latent laundering path; it does not paper over a live defect.

  1. DECEPTIVE TEST NAMES (3/3 fixed).
     (a) test_benchmarks.py::test_benchmark_correlation ->
         `test_benchmark_is_executable_and_every_measured_compound_is_detected`. It asserts
         no correlation coefficient, no ratio bound, no error bound — it is an executability
         + species-DETECTION smoke check, and detection is a coverage property, not an
         agreement one. Also TIGHTENED: two `pytest.xfail()` calls (imperative — they abort
         at that line) excused unexecutable benchmarks and incomplete coverage, i.e. exactly
         the properties under test. Measured: all 14 supported, all coverage 1.0, so both
         branches were dead; now assertions.
     (b) test_matrix_only_benchmark.py::..._is_quantitatively_supported ->
         `test_trikusuma_heated_pea_matrix_fit_is_recovered_to_within_1_05x`, with the
         tolerance explained: <=1.05x is expected BECAUSE the constants were back-fitted to
         this benchmark. Kept as a regression check on the fitting machinery; documented as
         NOT evidence of predictive accuracy.
     (c) THE POUVREAU CITATION IS FABRICATED. `test_pouvreau_pea_ph_family_...` cited
         "Pouvreau 2021, 55-65% higher headspace release at pH 4.5 vs 6.5". Verified against
         CrossRef (author search 2018-2024) and this repo's own
         citation_verification_ledger.md, where "Pouvreau" appears ZERO times: NO SUCH PAPER
         EXISTS. Three stacked defects: (i) ghost attribution — the real source, already
         repaired in data/lit on 2026-08-26 under a doi_repair, is Fischer, Cachon & Cayot,
         Food Res. Int. 2021, 150, 110760 (10.1016/j.foodres.2021.110760); (ii)
         over-generalisation — Fischer reports 59% for HEXANAL ONLY and explicitly "for
         example", so 2 of the 3 compounds asserted had no quantitative backing; (iii) wrong
         physical quantity — Fischer varied EXTRACTION pH (a manufacturing variable), the
         test varies pH as a release parameter at measurement time. And the 1.5-1.7 "band"
         was the model's own knob: src/headspace.py:241 gives exactly exp(2 x 0.235) =
         1.6000. Re-pinned as
         `test_ph_release_surrogate_reproduces_its_own_log_slope_knob`, a behavioural
         regression on an unanchored surrogate, citation claim REMOVED, `no_verifiable_source`
         note added, band TIGHTENED to the knob's exact value (a wide band around a
         deterministic closed form only hides drift). Also recorded: the three per-compound
         assertions were the same assertion three times (the factor is
         concentration-independent), and the pyrazine assertion is 1.0 by construction (a
         substring gate on "anal"/"enal"/"furan"), kept as a SCOPE guard and labelled as one.

  2. EMPTY-BODY TESTS: 15 removed, 1 implemented for real.
     tests/integration/test_cantera_sim.py — 13 `pass`-bodied tests across six classes,
     collecting and reporting PASS in 0.03 s. Not merely unimplemented but UNIMPLEMENTABLE:
     they invoke `--temp-profile linear:20to150:30min`, `--temperature`, `--time-total`,
     `--pH`, repeated `--output` for multiple formats, and import
     `src.cantera_integration.generate_mechanism_file` plus RMG fixtures. None exists — the
     module is src/cantera_export.py and the real CLI takes `--temp`/`--temp-ramp`/`--ph`.
     Uncommenting them would have tested a fictional program. Deleted with the reasoning
     recorded in-file, plus a note that the YAML-parses property they claimed IS covered for
     real by TestIsothermalSimulation. ONE was implementable and is now implemented:
     `test_cli_requires_the_precursors_argument` (subprocess, asserts non-zero exit and that
     stderr names the argument). File: 23 "tests" -> 9 real ones.
     tests/integration/test_temp_profiles.py::test_cli_ramp — wrote two fixture files, read
     neither, ended on `pass` under a comment admitting it tests nothing. Deleted; the ramp
     logic is genuinely covered by test_temperature_ramp above it.
     Residual `pass` sites checked and left: import guards, fake-class stubs, and
     except-clauses in parsing loops. All legitimate.

  3. SKIP / XFAIL AUDIT.
     KILLED — self-excusing (condition = "the thing under test is broken"):
       * test_blind_spots.py, both xfails. `pytest.xfail(...)` sat on the line BEFORE the
         assertion, so the assertion never ran and the test could never be retired by the
         code improving. `test_metal_catalysis_blind_spot` was worse: the xfail was its
         FINAL statement and the function contained no assertion at all. Both converted to
         declarative `@pytest.mark.xfail(strict=True)` with the assertions written out, so
         closing the gap now XPASSes -> FAILS -> forces promotion. Both gaps verified still
         real 2026-08-27 (matrix: both branches meaty 0.0; iron: bit-identical 0.11253319637981546).
         Recorded for the future implementer: the matrix comparison is DEGENERATE, not just
         negative — implementing inhibition alone will not flip it without a non-zero baseline.
       * test_error_handling.py x2 — `pytest.skip("Implementation may not validate barrier
         types")`. Also vacuous independently: step list `[]` means no barrier is ever looked
         up, so the test could not have detected validation either way. Both rewritten to
         assert measured behaviour (barriers: returns normally, NO boundary validation —
         stated as a known gap; concentrations: raises TypeError, but incidentally, from
         arithmetic rather than a validator — noted).
       * test_precursor_resolver.py x3 — `try: assert ...; except ValueError: pass`. Measured
         and pinned: thiamine/glutathione DO resolve; `resolve("glysine")` RAISES (there is no
         fuzzy matching — arguably the better contract, now asserted so adding one is
         deliberate); `resolve_many([])` returns `[]`.
     KILLED — dead (condition can never be true, file is tracked in git):
       test_uncertainty_propagation.py x2, test_fft_bottleneck.py, test_db_chemistry.py
       (module-level skipif), test_recommendation_engine.py (this one skipped the WHOLE
       pea/soy comparison, including the present file, if either was missing). All converted
       to assertions naming the tracked file.
     KEPT — legitimate, reason strings sharpened: tests/qm/* (xTB, MACE, PySCF, CREST, Sella
       binaries) and test_wave_i_tooling.py x2, where the precondition is NOT "someone forgot
       to run a generator" but that .gitignore excludes results/validation/* and
       family_sensitivity.json is not force-added. Consequence written into the file rather
       than hidden: the Wave I routing fix is NOT gated on CI. -> J2->J1 handoff.
     `xfail_strict = true` added to pytest.ini. Safe today: the suite contains exactly two
     xfail markers, both made strict explicitly above; there were no non-strict xfails left
     to surprise.

  4. NO-OP PARAMETRIZATION / ORPHANED SCAFFOLDING (tasks 2 + 6, same root).
     tests/benchmarks/ held ONLY `_lane_policy.py`, a helper for a Phase 3 QM authority lane
     whose tests were deleted 2026-04-21 (this ledger's finding S8). Nothing imported it; the
     directory collected zero tests. DELETED.
     That exposed a real no-op: test_skip_policy_registry.py spawned a pytest subprocess over
     tests/benchmarks and asserted skip_count == 0, lane_counts == {}, rows == [] — all true
     by construction, asserting that nothing is nothing, while the classification logic the
     module exists for was never touched. Replaced with a real test driving parse-and-classify
     through synthetic `pytest -rs` output (fast, deterministic, independent of which optional
     QM backends happen to be installed).
     SRC CHANGE (1 of 2, minimal, loudly logged): src/skip_policy_registry.py default lane
     list "tests/benchmarks", "tests/qm" -> "tests/qm". Every skip scan since April had been
     scanning an empty directory and reporting its zeros as a clean lane.
     tests/README.md updated: the removed lane is explained, and a "What counts as a
     legitimate skip" section now states the rule the audit applied.

  5. GUARD TESTS — tests/scientific/test_honest_headline_guards.py (NEW, 7 tests, 2.7 s).
     The durable defence: every honest headline lived only in English prose in README/AUDIT,
     which nothing enforces. Each guard RECOMPUTES a published number and cross-checks the
     doc mechanically, so docs and evidence cannot diverge silently. Pinned:
       * panel 14, strict-ready 0/14
       * evidence-role split 6/4/4 AND passes 0/6 predictive, 3/4 fit-recovery, 4/4 synthetic,
         7/14 aggregate. The SPLIT is pinned as well as the counts — otherwise reclassifying
         a benchmark from `predictive` to `fit_recovery` improves the headline while changing
         nothing about the model.
       * honest external-literature 1/3, not_evaluable 4, excluded_fitted_rows 2 AND
         excluded_fitted_rows_that_would_have_been_hits 2, median CI width 0.856 dex (pinned
         so coverage cannot be improved by widening intervals)
       * hold-out 0/5 genuine extrapolations PRE-WIDENING, matrix_sigma 2.0 asserted so the
         comparison stays the one the README claims, kind split 5/3, 8 matched points,
         median 32.79x, worst 2474.4x
       * pentose>>hexose 8.98x (ribose 981.3 / glucose 109.3), pinned as a two-sided band, NOT
         a floor — the existing `>= 3.0` floor test cannot detect the ratio climbing back
         toward the retired 15.8x
       * no_verifiable_source census, recomputed from scratch
     METHOD NOTE: most numbers are recomputed LIVE from tracked benchmark files rather than
     read from results/validation/*.json, because .gitignore excludes that directory wholesale
     — benchmark_summary.json and validation_overview.json DO NOT EXIST in a fresh clone. A
     guard reading them would skip itself into uselessness on CI, i.e. the exact pattern this
     wave removed. Full panel re-evaluates in ~1.3 s. prediction_uncertainty.json and
     external_validation_report.json are force-added and so are read from disk.

  6. FINDINGS.
     (a) THE POUVREAU CITATION IS FABRICATED — see 1(c). Not a test defect: a fabricated
         literature anchor that had been laundered into a "validation" test.
     (b) NO CHEMISTRY DEFECT behind the assert_balanced hole — measured, see 0.
     (c) The no_verifiable_source census was UNREPRODUCIBLE before this wave: no generator
         emits it, so README's numbers were unverifiable assertions. Now derived. The exact
         definition that reproduces the published triple: any object at any depth with
         source_status == "no_verifiable_source" over every parseable .json/.yml/.yaml under
         data/; "numeric" = the record contains a number anywhere. Yields 102 total / 80
         numeric / 62 in data/lit (= the runtime-reachable subset), matching README exactly.
         Independently corroborates J1's concurrent 84 -> 102 correction: 84 was data/lit only,
         and the 18-record difference is data/qm/phase33 (9) + phase35 (9). The 84/18 SPLIT is
         pinned in its own test, because a total of 102 can be preserved while data/qm quietly
         drops back out of scope and the honesty gain is undone.
     (d) src/skip_policy_registry.py was scanning an empty directory — see 4.

  7. J2->J1 HANDOFF (doc/config, not mine to edit).
     (a) GHOST CITATION STILL LIVE IN TWO PLACES. data/benchmarks/maillard_validation_benchmarks.md
         lines ~275-285 carry "Pouvreau 2021" as PRIMARY with the self-declaring placeholder
         DOI "(Approx. DOI: 10.1016/j.foodchem.2021.xxx — retrieve via Scopus)" and a fake data
         table (Hexanal ~340 vs ~205 ppb, "+59% at pH 4.5") — that table is where the test's
         205.0 input came from, i.e. invented absolute ppb reverse-engineered to give 59%.
         Line ~623 lists it in the "DOI still to retrieve" backlog. AND
         docs/reference/VALIDATION_CONTRACT.md lines ~165, ~192 still cite "the Pouvreau
         acidic-vs-less-acidic pea trend" as a live anchor. Correct both to Fischer, Cachon &
         Cayot 2021 / 10.1016/j.foodres.2021.110760, restrict the quantitative claim to
         hexanal, and note the extraction-pH vs release-pH mismatch. NOTE the citation gate
         cannot catch this class: there is no DOI to fail on, only a placeholder string.
     (b) results/validation/family_sensitivity.json is gitignored, so the two Wave I routing
         guards in test_wave_i_tooling.py never run on CI. Force-add it as was done for
         prediction_uncertainty.json / external_validation_report.json, and the skips become
         assertions. Same argument applies to benchmark_summary.json and
         validation_overview.json, which are the artifacts behind "14", "0/6" and "0/14" and
         are absent from every fresh clone — the guard tests work around this by recomputing,
         but the AUDIT/README tables cite artifacts a reader cannot open.
     (c) docs/reference/VALIDATION_CONTRACT.md line ~366 describes `tests/benchmarks/` as
         "intentionally skip-heavy today ... Phase 3 placeholders". That directory is now
         deleted (it had no tests at all, only an orphaned helper). Update or drop the line.

  8. RESIDUAL, NOT FIXED (honest list).
     * tests/unit/test_error_handling.py retains ~5 `try: assert ...; except: pass` blocks of
       the "either behaviour is acceptable" kind. Weaker than the cases fixed above — they do
       assert on the happy path — but they still absorb contract changes silently. Worth a
       pass; not done here.
     * The 0/6 guard pins `overall_status == "pass"`, matching how README/AUDIT count
       "without blocking gaps". Note that ALL 14 benchmarks carry non-empty `blocking_issues`;
       "pass" and "no blocking issues" are not the same predicate in this codebase. The
       stricter reading would be 0/14 everywhere.
     * The doc cross-checks assert token PRESENCE (e.g. "0/6 predictive" appears in AUDIT.md),
       not sentence structure, so they survive rewording but would not catch a number quoted
       correctly in one place and wrongly in another.

## Wave L1 — reaction-network topology vs. isotope-labelling literature (2026-08-27)

Deliverable: `docs/validation/isotope_topology_evidence.md` (new; the only file this wave
wrote besides this ledger entry). 28 model routes / route-absences scored against
isotope-labelling and pathway-elucidation studies: 8 CONFIRMED, 4 CONTRADICTED,
10 PARTIALLY SUPPORTED, 3 UNTESTED, 1 correct omission, 3 completeness gaps. Every citation
was retrieved during the review and tagged by access level (FULL / ABS / META / SECONDARY);
no verdict rests on a SECONDARY source, and 33 DOIs were re-verified against Crossref.

THE HEADLINE (CONTRADICTED, and it is the flagship). The primary MFT route installed by
Wave G1 — 1-deoxyosone -> NORFURANEOL -> +H2S -> MFT, cited to van den Ouweland & Peer
1975 — is contradicted by the decisive competition experiment. Cerny & Davidek 2003
(10.1021/jf026123f) spiked authentic norfuraneol into a [13C5]ribose/cysteine system: the
MFT came out "mainly 13C5-labeled, suggesting that it stems from ribose and that
4-hydroxy-5-methyl-3(2H)-furanone is UNIMPORTANT AS AN INTERMEDIATE", and they propose the
1,4-dideoxyosone route instead. Cerny & Davidek 2004 (10.1021/jf035265m) confirm
"1,4-dideoxypento-2,3-diulose as an intermediate of 2-methyl-3-furanthiol" positionally with
[1-13C]ribose. Independently, Hofmann & Schieberle 1998 (10.1021/jf9705983) — THE PAPER
`refit_sulfur_barriers_hofmann.py` FITS AGAINST, and this repo's only surviving sulfur
anchor — says in its own abstract that "Thiamin and norfuraneol/cysteine were less effective
precursors of MFT", with the highest MFT yield (1.4 mol %) coming from a C2+C3 recombination
(hydroxyacetaldehyde + mercapto-2-propanone) that the network cannot express at all.
MEASURED EXPOSURE: enumerating ribose+cysteine at pH 5.5 / 150 C, `Thiol_Addition_Norfuraneol`
is the SOLE producer of MFT — there is no second channel.
WHAT SURVIVES: the coarse topology is right. Cerny & Davidek show ribose fragmentation "did
not play a significant role" and the C5 skeleton stays intact, and the model's implicit atom
mapping (sugar C-1 -> 1-deoxyosone CH3 -> MFT 2-CH3) matches the observed
2-[13C]methyl-3-furanthiol from [1-13C]ribose. So this is a wrong-intermediate finding, not a
wrong-skeleton finding.
CANDIDATE FIX, and it REMOVES fiction rather than adding it:
    1-deoxy-2,3-pentodiulose C5H8O4 + 2[H] -> 1,4-dideoxy-2,3-pentodiulose C5H8O3 + H2O
    1,4-dideoxy-2,3-pentodiulose C5H8O3 + H2S -> MFT C5H6OS + 2 H2O    [EXACT, no 2[H]]
The sulfur-incorporation step needs NO reducing-equivalent token under the literature
topology; the reduction moves upstream onto a dehydration where an enediol donor is ordinary.
Norfuraneol should be KEPT as a species (Cerny & Davidek show it is the real precursor of
2-mercapto-3-pentanone) but demoted out of the MFT lane.
CONSEQUENCE FOR OPEN ITEM 1 (re-run the sulfur refit): do the route change FIRST. Refitting
`thiol_addition_norfuraneol` against Hofmann1998 before the route change fits the wrong route
more precisely, and the barrier is currently absorbing a route-selection error.

SECOND CONTRADICTION (lipid lane, live, cheap to fix). Nonanal is generated from the wrong
fatty acid. `predict_hexanal_generation` computes the whole hydroperoxide pool from
`linoleic_acid_pct` (lipid_oxidation.py:437) and emits 15 % of it as nonanal (:451); the
formulation path does the same at 12 % (:395). Miyazaki et al. 2023 (10.1093/bbb/zbac189,
read in FULL TEXT) resolves both linoleate hydroperoxide isomers: 13-HpODE -> hexanal, vinyl
hexanoate, 2-pentylfuran, 4,5-epoxy-2-decenal; 9-HpODE -> 2,4-decadienal, octanoic acid,
2-octenal, hexanal. Nonanal is in NEITHER list — it is the C9 fragment of the OLEATE double
bond. `LipidProfile.oleic_acid_pct` is declared, populated for pea and soy, carried through
the calibration registry, and READ BY NO CODE PATH IN THE REPO (grep-verified across src/,
scripts/, tests/); same for `alpha_linolenic_pct`. Hexanal + 2-pentylfuran co-production from
13-HpODE is supported, so the fix is local to the nonanal branch.

THIRD RISK (H4 residual, MEASURED). Wave I fix 8 gave MFT a cysteine-derived `[HH]` donor,
but the HEXOSE FURANONE step (reaction_templates.py:343, `elif can == hexose_do1 and h2 is
not None`) is still pool-gated on the token, and in a cysteine-free system the only producer
running before it is the pyrazine aromatisation. Measured on glucose+glycine, pH 5.5 / 150 C,
max_generations=3, by monkey-patching `_aminoketone_condensation` to return [] in-process
(no repo file touched): baseline 16 steps / 1 DMHF step / H2 producers
['Aminoketone_Condensation']; with the pyrazine lane off, 14 steps / 0 DMHF steps / no H2
producers. So PREDICTED FURANEOL FROM GLUCOSE IS STILL CONTINGENT ON PYRAZINE CHEMISTRY —
the exact failure mode Wave I fixed for MFT, surviving one lane over. It is worse than it
looks because Wang & Ho 2008 (10.1021/jf8012025) CONFIRM that route with a CAMOLA experiment
(1:1 [13C6]/[12C6]glucose -> 1:1 [13C6]/[12C6]DMHF = intact C6 skeleton; acetylformoin
detected as the precursor): the model is gating a confirmed route behind an unphysical
dependency. The cysteine/cystine donor is unavailable here, so this needs either the
sugar-derived enediol couple Wave I declined to invent, or the reduction lumped into the
cyclisation barrier with the token dropped.

OTHER CONTRADICTIONS. (a) HMF from FRUCTOSE goes through the 3-deoxyosone in this model
(Heyns -> hexose 3-deoxyosone -> HMF), and Perez Locas & Yaylayan 2008 (10.1021/jf8010245)
label-exclude it: "both fructose and sucrose showed much higher conversion rates than
3-deoxyglucosone thus precluding it as a major precursor of HMF in fructose and sucrose
solutions". Fructose is a registered precursor, so this is live. The glucose route is fine.
(b) ACRYLAMIDE's carbonyl-identity neutrality (already a repo open item, Wave I fix 18) is
contradicted quantitatively: Stadler et al. 2004 (10.1021/jf0495486) put the N-glycosyl of
Asn at ~2.4 mmol/mol vs 0.1-0.2 for the Amadori compound, and hydroxyacetone (>4 mmol/mol)
an order of magnitude above the alpha-dicarbonyls. The lumped `reducing_sugar_mM` cannot
express that.

WHAT HELD UP. The core sugar trunk survives isotope scrutiny intact — Yaylayan & Keyhani
2000/2001 (10.1021/jf000004n, 10.1021/jf000986w) independently reproduce every fragment the
model emits from it (pyruvaldehyde/acetol from intact C1-C2-C3 AND C4-C5-C6; glycolaldehyde
from C1-C2 and C5-C6), so `Retro_Aldol_Fragmentation` is CONFIRMED. FFT is the best-evidenced
route in the model: three independent results converge on furfural + H2S, including the
positional match ([1-13C]ribose -> 2-[13CH2]furfurylthiol, i.e. the sugar C-1 becomes the
CH2SH carbon exactly as the model routes it) and "furan-2-aldehyde/H2S showed a 10 times
higher efficiency". Acrylamide's Asn route is CONFIRMED by Zyzak et al. 2003
(10.1021/jf034180i) — and the two literature branches (loss of NH3 / loss of a substituted
imine) turn out to be split one each across the template layer and the curated layer.
DMHF/HEMF from pentose + amino acid is CONFIRMED with the exact C5+C1 / C5+C2 stoichiometry
the model uses (Blank & Fay 1996, 10.1021/jf950439o, with 13C-glycine and 13C-alanine). The
thiamine MFT cascade is CONFIRMED (Cerny & Briffod 2007, 10.1021/jf062874w), and its key
intermediate 5-hydroxy-3-mercapto-2-pentanone is a compound identified in exactly that system.
DMS's absence is a CORRECT omission (its precursor is S-methylmethionine, not methionine).

BONUS FOR THE PANEL. The repo's `thiamine_cys_glucose_120C_Bolton1994` benchmark (currently
pinned as a 724x honest failure) is sourced from Bolton, Reineccius, Liardon & Ba 1994,
"Role of Cysteine in the Formation of 2-Methyl-3-furanthiol in a Thiamine-Cysteine Model
System" — verified via Crossref to be a 34S-LABELLING STUDY. Secondary-source figures (NOT
retrieved; the chapter text is not online) say only ~7.5-8.0 % of the MFT carried 34S from
labelled cysteine, i.e. thiamine supplied nearly all the sulfur. IF that holds on inspection
of the chapter, the model's 724x miss is a magnitude problem in a lane whose TOPOLOGY IS
CORRECT — a materially different diagnosis from the sulfur-branch story the repo currently
tells, and checkable with one library visit. [P] flagged as the single highest-value
literature retrieval outstanding.

TOP 3 COMPLETENESS GAPS. (1) 2-ACETYL-1-PYRROLINE and the proline/ornithine lane: absent
(no precursors, no route), and the elucidating study is CAMOLA run INSIDE AN EXTRUDER at
135 C / 20 % moisture on a rice model — this repo's own process domain (Davidek et al. 2013,
10.1021/jf4004237: 2-AP by acylation of 1-pyrroline via C2 sugar fragments, major; C3, minor).
Both fragment donors are already species in the network. (2) LIPID-DERIVED STRECKER
INITIATION: 4-oxo-2-alkenals, 2,4-alkadienals and 4,5-epoxy-2-alkenals raise
phenylacetaldehyde yields 300-900 % with measured Ea 28-67 kJ/mol (Zamora & Hidalgo,
10.1021/jf301258s, 10.1021/jf305007y), and none of them are species; 2,4-decadienal, the
principal 9-HpODE product, is missing too. The model's only lipid x Maillard family makes
THIAZOLES instead. (3) The low-moisture C2+C3 recombination route to MFT (see the headline).
Runners-up: ~50 % of CML flux is the Amadori-oxidation route (Glomb & Monnier 1995,
10.1074/jbc.270.17.10017) and is unmodelled while CML sits 1204x low — worth trying before
any parameter is touched; 2,3-dimethylpyrazine and 2-ethyl-3,5-dimethylpyrazine are SCORED
TARGETS with no generative route (an isotope-traced route for EDMP exists,
10.1021/acs.jafc.9b07809); furfuryl alcohol and the odourless cysteine-S-conjugate sink
(Cerny & Guntz-Dubini 2013, 10.1016/j.foodchem.2013.04.043) are absent.

INCIDENTAL VERIFICATION. Wave I fix 4's flagship DOI repair is CORRECT: 10.1021/jf60199a045
resolves to Van den Ouweland & Peer (1975), "Components contributing to beef flavor. Volatile
compounds produced by the reaction of 4-hydroxy-5-methyl-3(2H)-furanone and its thio analog
with hydrogen sulfide", JAFC 23:501-505. Note the irony that the repaired-to-correct citation
now anchors a step the isotope literature retires.

FILES WRITTEN: docs/validation/isotope_topology_evidence.md and this entry. No src/, data/,
tests/, README or AUDIT file was touched; the two measurements were made in-process.

## Wave N — MFT route correction on isotope evidence (2026-08-27, orchestrator-executed)

THE FIX Pablo asked for ("can we find the right fix?"). Every step verified empirically
before landing; no barrier tuned; the honest numbers got worse and are reported as such.

WHAT CHANGED:
- `_norfuraneol_mft_route` -> `_furanone_and_mft_route` (src/reaction_templates.py).
  The norfuraneol + H2S + 2[H] -> MFT step is REMOVED; new steps:
    1-deoxy-2,3-pentodiulose + 2[H] -> 1,4-dideoxypento-2,3-diulose + H2O  [Deoxyosone_Reduction]
    1,4-dideoxypento-2,3-diulose + H2S -> MFT + 2 H2O   [Thiol_Addition_Pentodiulose, NO token]
  Both RDKit-verified exact. New species _PENTODIULOSE_14_DIDEOXY = CC(=O)C(=O)CCO
  (C5H8O3) in src/smirks_engine.py. Norfuraneol KEPT as a genuine furanone product
  (and documented as the demonstrated 2-mercapto-3-pentanone precursor, unimplemented).
- Same correction mirrored in the curated layer (src/curated_pathways.py pathway C).
- Barriers (src/barrier_constants.py): thiol_addition_norfuraneol 26.85 marked RETIRED
  (kept for provenance; the value was fitted THROUGH the contradicted route and does not
  transfer); thiol_addition_pentodiulose = 28.60 ESTIMATED UNCONSTRAINED (thiol_addition
  class value, deliberately NOT the 26.85); deoxyosone_reduction = 28.0 ESTIMATED
  (dehydration analogue, set equal to furanone_cyclisation to express NO prior preference
  at the 1-deoxyosone fork). Both have explicit FAST_BARRIERS keys because
  _canonical_fast_family would otherwise silently route Pentodiulose to the thiol_addition
  class and Deoxyosone_Reduction to the 45 kcal default (the known fallthrough bug class).
- ENGINE_FAMILY_LABELS extended (src/family_sensitivity.py); retired label kept (harmless
  superset, per test_engine_family_labels_cover_sources).

EVIDENCE (DOIs CrossRef-verified before writing; my first guess for the 2004 DOI was
WRONG — jf035125b — and caught by the check; correct is jf035265m):
- Cerny & Davidek 2003, 10.1021/jf026123f (JAFC 51:2714-21): spiking authentic
  norfuraneol into [13C5]ribose/cysteine leaves MFT mainly 13C5-labelled -> norfuraneol
  "unimportant as an intermediate"; 1,4-dideoxyosone route proposed.
- Cerny & Davidek 2004, 10.1021/jf035265m (JAFC 52:958-61): 1,4-dideoxypento-2,3-diulose
  positionally confirmed with [1-13C]ribose.
- Hofmann & Schieberle 1998, 10.1021/jf9705983: norfuraneol/cysteine among the LESS
  effective MFT precursors — the repo's own fit target contradicted its fitted route.
- Full dossier: docs/validation/isotope_topology_evidence.md (Wave L1).

MEASURED EFFECT (untuned): network enumeration confirms MFT sole producer is now
Thiol_Addition_Pentodiulose; norfuraneol has zero consumers. Hofmann1998: MFT 151.87 ppb
vs 342 = 2.25x under (was 1.45x — that number was carried by a barrier fitted through
the contradicted route); FFT 243.72 vs 200 = 1.22x over (was 1.10x). WORSE AND HONEST.
Guard tests fired exactly as designed (fold-band, ordering, headline guards) — being
reconciled by Wave M against regenerated artifacts, with dated causal comments.

TESTS RE-PINNED HERE (dated comments in each): test_chemistry_soundness (route tests ->
dideoxyosone route + retired-step-absent assertion), test_wave_h_2026_08 (RETIRED note
contract + new test_corrected_mft_route_barriers_are_estimates_not_fits pinning 28.60/
28.0 as ESTIMATED UNCONSTRAINED), test_wave_i_network_chemistry (family assertions).
168 chemistry/tooling tests green at handoff to Wave M.

[P] NEW OWNER DECISIONS:
1. Refit thiol_addition_pentodiulose against Hofmann1998? Now clean to do (single fit
   target, fit_target_gate will exclude that benchmark from scoring). Declined to do
   unilaterally: it re-couples the only sulfur anchor to a fitted constant.
2. Implement norfuraneol -> 2-mercapto-3-pentanone (the route Cerny & Davidek DID
   demonstrate for norfuraneol)? Adds a species; low cost, closes a sulfur sink.
3. The hexose legacy MFT shortcut remains non-mechanistic; the literature's dominant
   low-moisture route (C2 hydroxyacetaldehyde + C3 mercapto-2-propanone recombination,
   Hofmann's highest-yielding system, 1.4 mol%) is entirely absent from the network
   (mercapto-2-propanone is not a species). Wave L1 flags this as the likely dominant
   channel under extrusion conditions.

## Wave M — ONE regeneration + reconciliation pass over Wave N + Wave K/M (2026-08-27)

Scope: propagate two already-landed change sets. NOTHING in either was modified, and
nothing was tuned to recover a number that got worse.
  * CHANGE SET 1 (Wave N): MFT route corrected on isotope evidence.
    `Thiol_Addition_Norfuraneol` retired; `Deoxyosone_Reduction` (28.0 ESTIMATED) +
    `Thiol_Addition_Pentodiulose` (28.60 ESTIMATED UNCONSTRAINED) replace it.
  * CHANGE SET 2 (Wave K/M): six benchmark content corrections, each verified against the
    paper's own text and each carrying a `content_correction_note`.

### (a) REGENERATED HEADLINE NUMBERS — old -> new, with the one-line cause

  SNAPSHOTS. `refresh_internal_reproducibility_snapshots.py` NOTES bumped v6 -> v7 (it was
  already at v6 for Wave I, not v5) naming the Wave N route correction with both Cerny DOIs
  and the honest MFT/FFT movement. Run; the 4 pea/soy Internal2026 + ProtocolPilot2026
  snapshots re-pinned. MFT falls from rank 1 to rank 3 in the pea snapshot contract.

  HOFMANN1998 (cause: Wave N route correction)
    * MFT   235.32 -> 151.87 ppb vs 342   |  1.4533x under -> 2.2519x under
    * FFT   219.96 -> 243.72 ppb vs 200   |  1.0998x over  -> 1.2186x over
    * benchmark MALE 0.1019 -> 0.2192; its own contract (1.45x / 0.09 dex) UNTOUCHED.
    The 1.45x had been bought with `thiol_addition_norfuraneol` = 26.85, a barrier Wave H
    fitted THROUGH the contradicted route. Replacement ships at the un-fitted class value.

  PENTOSE >> HEXOSE (cause: Wave N; only the pentose side moved)
    * ribose MFT 981.3 -> 370.34 ppb; glucose 109.33 UNCHANGED; ratio 8.98x -> 3.3872x.
    * STRUCTURAL SHARE re-measured in-process (monkey-patched FAST_BARRIERS, no file
      touched): setting `thiol_addition_pentodiulose` (28.60) = `thiol_addition_hexose`
      (29.65) collapses the ratio to 1.1317x (was 1.29x under the old 26.85/29.65 gap).
      Also measured: additionally equalising `deoxyosone_reduction` gives 0.4926x, i.e. the
      pentose limb falls BELOW the hexose limb. So ~3x of the 3.39x rides on a 1.05 kcal/mol
      gap between two ESTIMATED barriers. Note 3.39x now sits close to the 3.0 floor in
      test_pentose_hexose_sulfur_ordering.py — one more correction of this size breaks the
      ordering claim outright, and that would be a finding.

  PANEL (cause: Wave K/M content corrections)
    * evidence-role split 6/4/4 UNCHANGED; passes: predictive 0/6 UNCHANGED,
      fit-recovery 3/4 -> 1/4, internal-synthetic 4/4 UNCHANGED, aggregate 7/14 -> 5/14.
    * status_counts pass 7 -> 5, scale-gap 7 -> 9.
    * inside-1.5x benchmark count 4 -> 1; outside 5 -> 8.
    * per-benchmark max_ratio:
        acrylamide_spi_extrusion_130C_ACSRef3      6.4236x -> 15.3871x (ref 62.62 -> 150)
        cys_ribose_140C_Hofmann1998                1.4533x ->  2.2519x (Wave N)
        pea_isolate_40C_PratapSingh2021            1.0025x ->  4.3661x (ref 260 -> 1138)
        soy_isolate_40C_PratapSingh2021            1.0007x ->  4.2690x (ref 380 -> 1621.71)
        resconi_2023_pbma_beef_identity_benchmark  3.1364x ->  4.5607x (ref 1040 -> 715.22)
      every other row bit-identical.
    * pea_isolate_40C_PratapSingh2021 ranking contract pass -> ORDER_MISMATCH: with the
      paper's real values hexanal (1138) outranks 2-pentylfuran (638) and the model has it
      the other way round (260.6 vs 638.3). A real ranking failure that was INVISIBLE while
      the reference agreed with the constant fitted to it.
    * validation_overview quantitative_point_count 41 -> 39 (the two hexanol rows).
    * strict-ready 0/14 and panel size 14 UNCHANGED.

  MC PANEL (prediction_uncertainty, n=200, seed 0)
    * benchmark_count 11, matched rows 35, aggregate coverage 29/35 — ALL UNCHANGED.
    * honest_literature_coverage 1/3, not_evaluable 4, excluded_fitted_rows 2/2, median CI
      width 0.8558 dex — ALL UNCHANGED.
    * only interval widths moved: fitted_row median CI 2.2752 -> 2.3205 dex,
      internal_synthetic 3.6499 -> 3.7657 dex (Wave N widened the MFT lane's spread).

  HOLD-OUT (cause: Wave K/M Li 2026 wrong-row correction; predictions byte-identical)
    * median fold error 32.79x -> 15.3074x; median |log10| 1.5158 -> 1.1849 dex.
    * pre-widening (ln-sigma 2.0): 3/8 -> 4/8; genuine extrapolations 0/5 -> 1/5.
    * shipped-sigma coverage 5/8 UNCHANGED; max_fold_error 2474.3839x UNCHANGED.
    * per point: 2-pentylfuran 49.835x -> 1.962x (meas 221.5 -> 5625.8, the Maltol row);
      nonanal 673.317x -> 272.626x (meas 29.42 -> 72.66, the Decanal row). The other six
      points are identical to 4 decimal places.
    * HOLD-OUT EXCLUSION RE-VERIFIED after regeneration: holdout_guard PASS (3/3 invariants);
      refit_projection_constants' own assertion that no external_validation path reaches the
      fit-target selector held; fit_target_gate PASS.

  OTHER ARTIFACTS regenerated with no headline movement: benchmark_index,
  experiment_value_ranking (top-5 order unchanged; rank-5 VoI 3.77 -> 3.87),
  loo_leverage, gap_heatmap (16 x 12, max VoI 13.36, 6 outside-CI cells),
  experiment_brief_cards, experiment_requests (--top 5, same five),
  validation_overview.{md,json,png} + family_validation_overview,
  family_implementation_status, hexanal_nonanal_{calibration_closure,resolution},
  objective_progress (+ its three PNGs). All figure generators ran WITH figures on the
  documented conda path (Wave I fix 16's mathtext fallback; dvipng still absent).
  `docs/slr_benchmark_evaluation.md` has NO generator and consumes none of the changed
  payloads — not regenerated, and its 41/37/27 numbers are unaffected.

### (c) THE REFIT MOVEMENT AFTER THE 4.3x ANCHOR CORRECTION — stated prominently

  MATRIX SIGMA (`derive_matrix_sigma_from_residuals.py`, leave-lane-out, hold-out never read)
    * n_residuals 6 -> 5 (the soy 1-Hexanol row is gone: the paper reports n.d.)
    * rms_ln_sigma        2.8620 -> 3.2520
    * bias_fold           3.4578 -> 3.9065
    * centered_sd_ln      2.8252 -> 3.3013
    * 90% CI            [1.9756, 5.4825] -> [2.1856, 6.7958]
    * shipped 2.86 still INSIDE the CI and was NOT moved.
    * the driver: soy ambient hexanal residual 2.2104x under -> 9.4283x under, because the
      measurement rose 4.27x while the prediction (171.91 ppb uncalibrated) did not.

  PROJECTION CONSTANTS (`refit_projection_constants.py`)
    * fitted reference_conversion_time_min  7943.28 -> 5011.87 min
    * objective at fit                      0.7327 -> 0.9408 dex
    * objective at SHIPPED tau (12589)      0.7396 -> 0.9599 dex
    * scored/total rows                     18/18 -> 16/16 (hexanol rows removed)
    * projection-sensitive-rows objective   0.7935 -> 0.8367 dex (5 rows)
    * FIT STILL NOT APPLIED. Unchanged decision, re-derived reasoning.
    THE OLD FIT WAS RECOVERING WRONG DATA. This is the load-bearing sentence of the wave.
    At the SHIPPED constants the pea/soy ambient hexanal rows read 260.65 and 379.88 ppb —
    i.e. the model reproduces 260 and 380, the numbers the observability factors were
    back-solved from, to 4 significant figures. Those numbers are not in the paper. So the
    0.7396 dex objective of record was ~0.00 dex on two rows that were agreeing with a
    transcription error; the honest values put them at 0.6401 and 0.6303 dex.
    GENERATOR BUG FIXED WHILE HERE (no numeric effect on the fit): `wave_h_decision.
    why_not_applied` HARD-CODED the two residuals the decision turns on ("MFT collapses to
    0.048 dex (from 0.75)", "furfural grows to 1.281 dex, i.e. ~19x OVER", "~1050 ppb",
    "78% of the pool", "5.6x allocation gap", "5250 ppb"). Every one was stale. They are now
    read out of the rows the run itself computes. Current, derived: MFT 0.3526 -> 0.0468 dex
    at the optimum; Resconi furfural 0.6590 -> 1.0574 dex, i.e. 4.6x -> 11.4x OVER; budget
    scale 2.51x. The allocation arithmetic was re-measured 2026-08-27: the Hofmann system's
    five tracked species total 1126.37 ppb against a measured MFT+FFT of 542 ppb, furfural
    taking ~51% (not ~78%). The record also now says out loud that
    `sulfur_barrier_refit_hofmann.{json,md}` is STALE — it profiles a retired family.

### (b) EVERY RE-PINNED TEST AND ITS DATED CAUSE

  18 failures at entry, 0 at exit. Every one was a stale pin; no failure was "fixed" by
  touching src/ or data/ except the two propagation items called out at the end.

  1. tests/scientific/test_free_aa_quantitative_regression.py
     MFT band (1.20, 1.453, 1.90) -> (1.85, 2.252, 2.80); FFT (1.00, 1.100, 1.45) ->
     (1.00, 1.219, 1.60). Cause comment: Wave N route correction, both Cerny DOIs, plus an
     explicit "the numbers got worse and were NOT clawed back" paragraph. Header block's
     stale "max_ratio 1.4533 / MALE 0.1019" updated to 2.2519 / 0.2192. Bands NOT widened in
     relative terms.
  2. tests/scientific/test_honest_headline_guards.py x3
     (a) fit-recovery passes 3 -> 1, aggregate 7/14 -> 5/14, docstring rewritten. Cause:
         Pratap-Singh content correction. Failure message now names the single survivor
         (Trikusuma2019) and warns its source is still content-unverified.
     (b) hold-out guard renamed `..._1_of_5_...`; genuine_extrapolation_hits 0 -> 1;
         median_accuracy_fold 32.79 -> 15.31; README token "0/5" -> "1/5". max_fold_error
         2474.4 deliberately UNCHANGED (the extreme misses did not move) — that is the
         evidence it was a reference correction, not a model change.
     (c) pentose guard renamed `..._is_3_39x_not_the_retired_8_98x_or_15_8x`; ribose
         981.3 -> 370.3, ratio 8.98 -> 3.39, doc token "8.98" -> "3.39". Glucose assertion
         kept at 109.3 with a note that it was NOT touched by Wave N, which is what
         identifies the cause.
  3. tests/scientific/test_deep_research_benchmark_wiring.py x2
     acrylamide band [4.5, 9.0] around 6.42x -> [11.0, 20.0] around 15.39x; Resconi furfural
     [2.5, 4.0] around 3.186x -> [3.6, 5.7] around 4.561x. Both CONTENT corrections, both
     with the prediction unchanged (9.7484 ppb and 3261.91 ppb). Band widths held at the
     same RELATIVE span as before, so neither was loosened.
  4. tests/scientific/test_benchmark_summary.py x2
     pea ranking_contract_status "pass" -> "order_mismatch". And
     `test_..._markdown_includes_gap_labels` asserted the "matrix-only intake path is
     executable" note, which the summary only emits when a row has NO other blocking issue;
     both Pratap-Singh rows now carry real scale gaps, so it asserts the labels they
     actually produce ("max ratio", "matrix ranking contract: order_mismatch") — which is
     what the test's own name claims.
  5. tests/scientific/test_matrix_assertions.py
     pea overall_status/adverse_order_status "pass" -> "fail", plus new explicit pins on
     ranking_contract_status == order_mismatch, max_ratio ~4.366 and ratio_status == fail.
     TIGHTENED, not relaxed.
  6. tests/scientific/test_matrix_headspace_ph_validation.py
     The blanket `ratio <= 1.01` over every compound of both benchmarks is SPLIT: the
     2-pentylfuran rows keep <=1.01 (their measured values were verified verbatim, so the
     recovery is still real), hexanal is pinned two-sided at 4.366x / 4.269x with a
     direction assertion. Two-sided is STRICTER than the old one-sided bound.
  7. tests/scientific/test_matrix_only_benchmark.py x5
     comparisons 3 -> 2 (hexanol removed: paper says n.d.); the three-way ordering assertion
     rewritten to assert the measured ordering AND the predicted ordering separately, so the
     pea DISAGREEMENT is asserted rather than tolerated; hexanal limits 1.25/1.10 (pass
     bands over an algebraic recovery) -> two-sided 4.366/4.269 pins; soy contract targets
     and adverse_markers lose "hexanol"; the `fitted_to_benchmark` comment extended to say
     the values those factors were solved from are not in the paper.
  8. tests/unit/test_audit_remediation_carried_2026_08.py x2 — see PROPAGATION below.

  NOT RE-PINNED, DELIBERATELY: `cys_ribose_140C_Hofmann1998`'s own validation contract
  (1.45x / 0.09 dex), the acrylamide contract (1.5x / 0.2 dex), the Resconi tolerance, the
  Pratap-Singh ratio tolerance (2.0x) and ranking contract, the 3.0x pentose/hexose floor.
  Every one of these now FAILS visibly in the regenerated artifacts. That is the point.

### PROPAGATION EDITS (not test re-pins; each changes zero numbers)

  * scripts/generators/generate_family_implementation_status.py — its unmapped-family
    ASSERTION fired correctly on `Deoxyosone_Reduction` and `Thiol_Addition_Pentodiulose`.
    Both added to `_FAMILY_TO_SLR` under `amino_acid_sugar_core` (the same limb as the step
    they replace); the retired key kept so an older tree still classifies.
  * scripts/generate_reaction_network.py — `uncoloured_reaction_families()` and
    `stale_colour_map_keys()` both fired. `Thiol_Addition_Pentodiulose` INHERITS the retired
    step's #20B2AA so the MFT lane keeps its hue across the route change;
    `Deoxyosone_Reduction` gets #6B8E23, adjacent to `Furanone_Cyclisation` because the two
    are the competing branches of the same 1-deoxyosone fork. Retirement recorded in the
    file's own "retired keys" comment block. NOTE: `docs/assets/reaction_network.pdf` was
    NOT regenerated (untracked, and outside this wave's generator list).
  * data/protocols/example_matrix_experiment_intake.yaml — a DOCUMENTATION DRIFT FIXTURE
    whose own notes say its ranks track the model and must be refreshed after each
    intentional model change (Wave I refreshed it for the same reason). Wave N drops MFT
    below both FFT and 2,5-dimethylpyrazine, so ranks 3/4/5 are now FFT / DMP / MFT. Without
    this the support-delta blocker reads "ranking contract not yet passing", which would be
    a fixture artefact reported as a promotion blocker.
  * src/matrix_calibration_registry.py — COMMENT ONLY, no constant touched. Its Wave I block
    asserted "the anchor is a real measurement from a real paper", which is no longer true
    for the two hexanal rows and never was for the hexanol rows. Now carries the paper's
    real values, the reason the factors were NOT refitted, and the consequence (the lanes
    under-predict by exactly the correction factor).

### (d) GENUINE NEW BUGS FOUND — reported, not fixed

  1. THE MATRIX OBSERVABILITY FACTORS ARE FITTED TO VALUES THAT ARE NOT IN THE PAPER.
     `MATRIX_BENCHMARK_BASE_MARKER_YIELDS` + the pea/soy ambient registry entries were
     back-solved from hexanal 260 / 380 and hexanol 80 / 120. Molecules 2021, 26, 4104
     Table 1 gives 1138.00 / 1621.71 and n.d. / n.d. The 2-pentylfuran factors are fine
     (638 / 2492 verified verbatim). Consequence: the ambient lane under-predicts hexanal by
     exactly 4.37x / 4.27x, the hexanol entries anchor to a compound the paper says was not
     detected, and the pea ranking contract now genuinely fails. REFITTING IS AN OWNER
     SCIENCE DECISION and was not done — doing it in the same pass as the Wave N chemistry
     change would make the resulting agreement unattributable. [P]
  2. `results/validation/sulfur_barrier_refit_hofmann.{json,md}` IS STALE, NOT WRONG-BY-DRIFT:
     it profiles `thiol_addition_norfuraneol`, a family no step emits any more, and README /
     AUDIT / test comments cited it as a live finding. Flagged in the projection record and
     in AUDIT.md rather than re-run: re-running it would fit the one remaining sulfur anchor,
     whose own mol%->ppb conversion is undocumented (Wave K finding 6). [P]
  3. `refit_projection_constants.py` hard-coded the residuals its DECISION PROSE quotes, so
     the artifact kept publishing Wave H numbers through two chemistry changes and a data
     correction. FIXED HERE (derived from the run's own rows). The same pattern should be
     assumed elsewhere: any generator that writes a narrative alongside numbers is a
     candidate. [P]
  4. README.md stated "Only 3 of the 8 [hold-out] points test extrapolation at all. The
     other 5 ... are reproducibility comparisons." The artifact says the opposite and always
     has (5 genuine_extrapolation, 3 in_panel_rescoring). Pre-existing, not caused by this
     wave; corrected in the doc sync.
  5. README's parity figure caption said "across 16 benchmarks" against a 14-benchmark panel;
     its "coverage fell to 18%" was the retired 2/11 accounting; "26 of the 41 matched rows"
     mixed the MC panel's 26 with validation_overview's 41. All three corrected.

### DOC SYNC

  README.md: fitted/synthetic median CI widths 2.28/3.65 -> 2.32/3.77; fit-recovery 3/4 ->
  1/4 and 7/14 -> 5/14 with a new pull-quote on WHY a fit-recovery that stops recovering is
  the cleanest possible proof it was never evidence; a new paragraph on the Wave N route
  correction with both Cerny DOIs and the honest 1.45x -> 2.25x; projection objective
  0.74 -> 0.96 dex with the reason (the references got more accurate, the model did not
  change); hold-out row 0/5 -> 1/5 and 32.79x -> 15.31x with a new block naming the wrong
  table rows; the 3-of-8 / 5-of-8 inversion fixed; pentose block 8.98x -> 3.39x with the
  structural share 1.29x -> 1.13x and the "what improved is the KIND of unconstrained"
  point; family-08 lane 6.4x -> 15.4x; parity caption 16 -> 14; and a new block
  "The matrix over-prediction story, re-derived on corrected anchors" that states plainly
  that the over-prediction did NOT largely disappear (median halved, one 49.8x point became
  1.96x, the 2474x/1117x/273x extremes untouched) while a NEW 4.3x under-prediction appeared
  in-panel where there had been perfect agreement.
  AUDIT.md: table rows for blocking gaps and hold-out; the pentose ordering sentence; a
  Round-3 addendum on the superseded Hofmann numbers; the "real norfuraneol MFT route ...
  selectivity is now structural" bullet corrected (it was wrong when written); Open item 1
  rewritten (the orphaned constant was resolved by retiring the step, and the refit record is
  now stale); and a new section "Round 3: content verification + route correction" with the
  6 findings table, the ~25% base rate, Europe PMC fullTextXML as the access route that beat
  the MDPI 403s, the isotope evidence, a Round-3 cost table, and the generalisation that the
  live panel is the BEST-audited part of the repo, not the worst.
  docs/reference/VALIDATION_CONTRACT.md §3E hold-out numbers 0/5 -> 1/5, 32.79x -> 15.31x,
  3/8 -> 4/8, with the wrong-row cause spelled out.

### (e) GATES + SUITE

  citation_gate PASS — 79 files, 904 DOI-bearing fields, 310 unique DOIs (was 74/846/288
  at Wave I). WAIVERS and TEXT_SURFACE_WAIVERS are both still the empty tuple. The new
  Cerny DOIs (10.1021/jf026123f, 10.1021/jf035265m) and the PMC identifiers passed without
  any formatting change being needed.
  holdout_guard PASS — 3/3 invariants, re-run after regeneration.
  fit_target_gate PASS — both checks.
  FINAL SUITE (tests/unit + tests/scientific + tests/integration + tests/scripts, documented
  conda path): **1242 passed, 1 skipped, 2 xfailed, 0 FAILED** in 809.72 s. Entry state was
  1224 passed / 18 failed. The 1 skip and both xfails are the declared, strict-marked ones
  from Wave J2 (xfail_strict = true is on, so neither xfail can silently start passing).
  All figure generators ran WITH figures; nothing was run with --skip-figures.

  NOT COMMITTED, NOT STASHED — handed to the orchestrator as instructed.

### [P] CARRIED FORWARD FROM THIS WAVE

  1. Refit the pea/soy ambient matrix observability factors against the corrected
     Pratap-Singh anchors (owner science decision; the constants currently reproduce values
     that are not in the paper). Note the fit_target_gate implication: they are already
     `fitted_to_benchmark`, so refitting changes the size of the recovery, not its status.
  2. Re-run `refit_sulfur_barriers_hofmann.py` against the corrected route — blocked behind
     recovering the Hofmann full text and writing down its mol%->ppb conversion.
  3. Audit every generator that writes narrative prose next to computed numbers for the
     hard-coded-residual pattern fixed in `refit_projection_constants.py`.
  4. `docs/assets/reaction_network.pdf` needs regenerating for the new family colours.
  5. ~200 content-unverified anchors in `data/lit/` sit behind a ~25% measured content-error
     base rate. Trikusuma 2020 is the highest-value single retrieval: it is the last
     fit-recovery benchmark still scoring `pass` and the last unchecked pillar of the matrix
     lane.

## Commit assembly note (2026-08-27, orchestrator)

- data/lit/timeseries/ (4 YAML files, 712 points, Wave L2 partial output) is DELIBERATELY
  NOT COMMITTED: the harvesting agent was killed before its verification report; the
  files carry honest access labels (figure_readoff/full_text/mixed) but nobody has
  checked a single point against the papers, and one file lacks a DOI. Policy: no
  unverified numbers enter the repo. [P] a verification pass over these four files
  (Martins 2003/2005, Brands) unlocks the first real kinetic calibration surface.
- Wave J1b/J2b (killed by a user interrupt mid-wave) turned out to have LANDED most of
  their scope before dying, and Wave M's certifying suite (1242/0) ran with it all in
  tree: .gitignore un-ignores WITH measured before-number (5 clean-checkout failures,
  matching CI's own historical record), data/qm + ingest_templates + Gemini_Deep_Research
  staged (68 files), Gemini warning README present, ci.yml strict-benchmarks advisory job
  + honest not-covered block, pytest.ini xfail_strict with both xfails made strict,
  smirks balance-test hardening. NOT done from their briefs: post-un-ignore clean-checkout
  AFTER-verification (running now, pre-push gate), citation-gate extension to data/qm,
  deceptive-test renames beyond the smirks files, empty-body test cleanup audit,
  tests/benchmarks scaffolding cleanup. Carried as [P] wave J follow-up.

## Wave O — refit the matrix observability factors onto the CONTENT-CORRECTED anchors (2026-08-27)

Owner-approved. Closes [P] item 1 carried forward from Wave M. Scope: the observability
factors and nothing else — no barrier, no projection constant, no marker yield.

### (a) THE HOLD-OUT, OLD -> NEW — the headline of this wave, and it is a REGRESSION

  `results/validation/external_validation_report.json`, n=200 seed 0, hold-out never fitted.

  | point                                   | measured | p50 old  | p50 new  | fold old  | fold new  |
  | --------------------------------------- | -------- | -------- | -------- | --------- | --------- |
  | bi_2020_raw_pea / hexanal               |   1260.0 |    234.66 |   1013.0 |    5.3695 |    1.2437 |
  | bi_2020_roasted_pea / hexanal           |    324.0 |  801700.4 | 801700.4 | 2474.3839 | 2474.3839 |
  | li_2026_hme / 1-hexanol                 |    20.04 |   22394.3 |  22394.3 | 1117.4799 | 1117.4799 |
  | li_2026_hme / 2-pentylfuran             |   5625.8 |   11038.6 |  11038.6 |    1.9621 |    1.9621 |
  | li_2026_hme / hexanal                   |    605.6 |   13066.2 |  56409.0 |   21.5756 |   93.1471 |
  | li_2026_hme / nonanal                   |    72.66 |   19809.0 |  19809.0 |  272.6257 |  272.6257 |
  | liu_2023_ppi / hexanal                  |    51.96 |    234.66 |   1013.0 |    4.5161 |   19.4972 |
  | liu_2023_ppi / nonanal                  |    15.81 |    171.70 |   171.70 |   10.8603 |   10.8603 |

  * median |log10| 1.1849 -> 1.6296 dex; median fold 15.3074x -> **42.6159x**
  * shipped-sigma coverage 5/8 -> **4/8** (0.625 -> 0.500)
  * pre-widening (ln-sigma 2.0) coverage 4/8 -> 3/8
  * pre-widening genuine_extrapolation_hits **1/5, UNCHANGED**; max_fold_error **2474.3839x,
    UNCHANGED**. Shipped-sigma genuine-extrapolation hits 2/5 -> 1/5.
  * Three points moved. Five are byte-identical.

  WHY A REFIT TO VERIFIED VALUES MAKES THE HOLD-OUT WORSE. The pea ambient lane carries two
  external measurements of nominally the same system that disagree by 24x: Bi 2020 at
  1260 ppb and Liu 2023's band midpoint at 51.96 ppb (band 15-180). The erroneous 260 ppb
  the old constants reproduced sat almost exactly at their geometric mean —
  sqrt(1260 x 51.96) = 255.9 against a shipped prediction of 234.66. Being wrong in the
  middle of a contradiction scores better than being right at one end of it. The verified
  anchor (1138 ppb) agrees with Bi to 1.11x and sits 6.3x above the top of Liu's band. Which
  is representative of commercial PPI is a literature question this repo cannot settle; note
  Liu's number is a geometric midpoint WE constructed and its band was never covered under
  the old constants either (234.66 was already 1.3x above band max).
  The Li 2026 movement is NOT a fit: it comes through the soy heated_matrix hexanal factor,
  which is DEFINED as `soy ambient baseline x (1 - 0.7060)` (Shu 2024 attenuation). Freezing
  it would have left a corrected fit composed with a stale baseline — the exact defect this
  wave exists to remove — so it was propagated.

### (b) OLD -> NEW FACTORS, WITH PROVENANCE

  Generator: `scripts/generators/refit_matrix_observability_pratap_singh.py`
  Record: `results/validation/matrix_observability_refit_pratap_singh.{json,md}` (both
  force-added in .gitignore; the registry cites them from the constants and
  `fit_target_gate` reads them).

  | lane                              | role       | before      | after     | move     |
  | --------------------------------- | ---------- | ----------- | --------- | -------- |
  | pea_iso/ambient_slurry/hexanal    | fitted     | 1.0         | 4.31725   | 4.3172x  |
  | soy_iso/ambient_slurry/hexanal    | fitted     | 0.453/0.205 | 9.54007   | 4.3172x  |
  |                                   |            | = 2.2097561 |           |          |
  | soy_iso/heated_matrix/hexanal     | propagated | 0.6496683   | 2.80478   | 4.3172x  |

  Each carries `previous_value=`, a dated causal comment and a pointer to the record.
  Fitted against: pea 1138.00 +/- 297.30 ppb and soy 1621.71 +/- 159.69 ppb, Molecules 2021,
  26, 4104 Table 1 via Europe PMC PMC8271896 (Wave K).

### (c) FIT QUALITY, AND WHAT DID *NOT* SATURATE

  ONE free parameter — a shared multiplicative scale on BOTH ambient hexanal factors — not
  two. Two factors against two rows is exactly determined: the residual would be 0 by
  arithmetic and would measure nothing. The one-parameter fit is reported alongside the
  two-parameter alternative (`alternative_per_lane_fit`, objective exactly 0.0, declined).

  * objective (sum of squared log10 residuals over the 2 anchored rows):
    **0.807024 -> 0.000048 dex^2**
  * adopted shared scale **4.317249x**; pea alone wanted 4.36606x, soy alone 4.26899x
  * residual after the fit: **1.0113x on BOTH rows** — pea reads just UNDER its anchor,
    soy just OVER, because a shared scale must bracket the two required scales
  * MUTUAL CONSISTENCY: the two corrected anchors agree to 1.1%. The transcription error was
    a common ABSOLUTE-SCALE error; the pea-vs-soy release ratio the registry encodes
    (2.2098) survives it untouched (a per-lane fit would have moved it to 2.1606).
  * NOTHING SATURATED. Search range [1e-3, 3e1] (4.5 decades, deliberately wider than the
    full span of the registry's own factors, 0.0095957-5.9203). The optimum sits 3.64 dex
    above the floor and 0.84 dex below the ceiling. `bounds_check.hit_a_bound = false`.
  * LINEARITY MEASURED, NOT ASSUMED: predictions re-evaluated at s and 2s give a ratio of 2
    to 1e-6, so J(s) is an exact parabola in log10(s) and the closed form is legitimate.
  * REPORTED, NOT FITTED (implied scales all within 5e-4 of 1.0, i.e. below the 0.01 dex
    materiality threshold): pea/soy ambient 2-pentylfuran (638 / 2492 verified VERBATIM —
    moving them would be fitting rounding) and all three Trikusuma heated-pea factors
    (that benchmark was not content-corrected and is still the last unread pillar of the
    matrix lane).
  * UNANCHORED, so unfittable: the 1-hexanol factors. The paper reports n.d. in both
    matrices, so the shipped soy value `0.143/0.063` is a ratio of two numbers (120 and
    80 ppb) that appear NOWHERE in it. Left untouched and flagged — it is the constant
    behind the hold-out's 1117x 1-hexanol miss. [P]
  * The pea/soy CLASS-level anchors (aldehyde 1.0 / 2.209) were deliberately NOT moved: they
    are the fallback for every aldehyde in every lane and no benchmark constrains that
    transfer. Because the refit preserved the soy/pea ratio, 2.209 remains correct AS A
    RATIO. [P]

### (d) MATRIX SIGMA RE-DERIVATION — VERDICT: 2.86 STILL DEFENSIBLE, AND THE REFIT SAYS NOTHING

  `derive_matrix_sigma_from_residuals.py` re-run after the refit produced a **BIT-IDENTICAL**
  artifact: n=5, rms_ln_sigma 3.2520, bias_fold 3.9065, centered_sd 3.3013, 90% CI
  [2.1856, 6.7958], shipped 2.86 INSIDE. Not luck — structural: `_uncalibrated_prediction_ppb`
  multiplies the oxidation load by the base MARKER YIELD and never reads an observability
  factor, because the uncalibrated tier exists to describe a lane that has NO such
  calibration. Consequence written into the generator's header: **no refit of these constants
  can ever be used to justify narrowing this prior.** The soy ambient hexanal residual stays
  9.4283x under in the uncalibrated view even though the calibrated view now reads 1.0113x.

### (e) OTHER REGENERATED NUMBERS

  BENCHMARK PANEL (`benchmark_summary.json`)
    * pea max_ratio 4.3661 -> 1.0113; MALE 0.3201 -> 0.0025; scale_status fail -> pass;
      overall_status scale-gap -> **pass-no-ranking**; ranking_contract_status
      order_mismatch -> **pass** (the model now puts hexanal above 2-pentylfuran, matching
      the paper — as FIT RECOVERY: the winning compound is the one whose factor was solved
      from its own measurement).
    * soy max_ratio 4.2690 -> 1.0113; MALE 0.3152 -> 0.0025; same status movement.
    * status_counts: scale-gap 9 -> 7, pass 5 -> 5, pass-no-ranking 0 -> 2.
    * evidence-role split 6/4/4 UNCHANGED; predictive passes **0/6 UNCHANGED**;
      fit-recovery passes **1/4 UNCHANGED**; aggregate **5/14 UNCHANGED** — because
      `pass-no-ranking` is not `pass`. Now pinned explicitly in the headline guard so the
      improvement cannot land invisibly.
    * strict-ready **0/14 UNCHANGED**; both rows still `strict_gate_blocked`.
  MC PANEL (`prediction_uncertainty.json`): UNCHANGED — the matrix-only benchmarks are not
    in it (11 of 14), so honest_literature_coverage 1/3, not_evaluable 4,
    excluded_fitted_rows 2/2 and median CI 0.8558 dex all hold.
  VALIDATION OVERVIEW: experimental benchmarks inside 1.5x **1 -> 3**, outside 1.5x 7 -> 5;
    worst quantitative point unchanged (CML 1203.68x). family_validation_overview SLR-02
    mean |log10| 0.085 -> 0.001.
  PROJECTION REFIT (`refit_projection_constants.py`, still NOT applied): objective at the
    shipped tau 0.9599 -> **0.8811 dex**, at the fit optimum 0.9408 -> **0.8620 dex**; fitted
    tau UNCHANGED at 5011.87 min. THE 0.079 dex GAIN IS ENTIRELY THE TWO FIT-RECOVERY ROWS
    (pea 0.6401 -> 0.0049 dex, soy 0.6303 -> 0.0049 dex) AND IS NOT EVIDENCE.
  SNAPSHOTS: `refresh_internal_reproducibility_snapshots.py` NOTES v7 -> v8. The four
    synthetic snapshots did NOT move at all — see the finding below.

### (f) GENUINE FINDING FOUND WHILE HERE (reported, not fixed)

  THE MATRIX OBSERVABILITY REGISTRY IS UNREACHABLE FROM THE RECOMMEND PATH. On the
  `matrix_precursor_augmented` lane the network species arrive labelled by SMILES
  (`CCCCCC=O`), not by name, so `describe_matrix_calibration` finds no record and returns
  `calibration_observable_factor=None`, applying 1.0. Measured directly on
  `soy_isolate_ribose_cysteine_100C_45min_Internal2026`: Hexanal metadata reads
  `calibration_factor: 1.0, calibration_evidence_strength: heuristic,
  calibration_fallback_mode: class_level`, and the snapshot is bit-identical after a 4.32x
  change to the soy heated hexanal factor. Two consequences: (i) the internal snapshots are
  insensitive to every constant in `src/matrix_calibration_registry.py`, so they cannot
  detect drift in it; (ii) the process_state is recomputed there as `intermediate_matrix`
  (100C/45min/aw0.95) rather than the file's declared `aqueous_pre_extrusion_model`. Recorded
  in the snapshot generator's NOTES. [P]

### (g) TESTS RE-PINNED (each with a dated causal comment)

  1. `test_matrix_headspace_ph_validation.py` — `_PRATAP_SINGH_EXPECTED_RATIOS` hexanal
     4.366/4.269 -> 1.0113/1.0113. The blanket "predicted < measured" direction assertion is
     SPLIT per lane (pea under, soy over) because a shared scale must bracket the two
     required scales — STRICTER than before, and it is the fingerprint that would break if
     someone re-added per-lane freedom. 2-pentylfuran's <=1.01 branch untouched.
  2. `test_matrix_only_benchmark.py` — same ratio re-pin; the ordering assertion rewritten
     from a hard-coded "2-pentylfuran first" to "predicted ordering == measured ordering",
     which now holds on BOTH lanes (fit recovery on the pea side); direction split per lane;
     the `fitted_to_benchmark` label assertion left UNCHANGED with a note that its stability
     across the refit is the point.
  3. `test_matrix_assertions.py` — pea ranking_contract order_mismatch -> pass, max_ratio
     4.366 -> 1.0113, ratio/adverse_order/overall fail -> pass. `strict_gate_blocked` stays
     True and is annotated: a fit-recovery row cannot promote itself.
  4. `test_benchmark_summary.py` — pea ranking_contract order_mismatch -> pass, plus NEW
     explicit `evidence_role == "fit_recovery"` assertions on both rows so the restored pass
     cannot be read as evidence. The markdown gap-label assertion moves off the
     Pratap-Singh rows (their gaps closed) onto Hofmann1998, which still has a real one.
  5. `test_honest_headline_guards.py` x2 — hold-out median 15.31 -> 42.62 with the full
     3-point movement table and the geometric-mean explanation, plus new pins on
     ci_coverage_hits 4 / rate 0.5; and the 0/6 guard gains a `pass-no-ranking` census
     (0 -> 2, both named, both asserted to be `fit_recovery`) so "the guard did not move"
     cannot hide a status improvement.
  6. NEW `tests/scientific/test_matrix_observability_refit_wave_o.py` (3 tests) — shipped
     constants vs the published refit to 5e-4 dex; the 1.0113x residual and the
     one-parameter declaration; and the record's `per_row_recovery` fit-target declaration.

  NOT RE-PINNED, DELIBERATELY: the Pratap-Singh ratio tolerance (2.0x), the strict-gate
  exclusion, the 0/6 / 5/14 / 0/14 headlines, `max_fold_error`, and the pre-widening 1/5.

### (h) DOC SYNC

  README.md: fit-recovery block gains a Wave O paragraph (shared scale, 1.0113x residual,
  status unchanged, ratio preserved); hold-out surface row 15.31x -> 42.62x and "2/5 via the
  wider prior" -> 1/5 at both priors; a new block on the 24x Bi-vs-Liu contradiction and the
  geometric-mean coincidence; the Wave M "matrix over-prediction re-derived" block updated
  with the Wave O median and the reason the sigma derivation is invariant.
  AUDIT.md: hold-out row; a new "Wave O — the refit, and the price it charged" section under
  Round 3 with the 8-point table, the one-vs-two-parameter argument, the unanchored hexanol
  finding and the sigma verdict; three new rows in the "What Round 3 cost" table.
  docs/reference/VALIDATION_CONTRACT.md §3E: 15.31x -> 42.62x, 4/8-5/8 -> 3/8-4/8, and the
  contradiction spelled out.

### (i) GATES + SUITE

  citation_gate PASS — 79 files, 904 DOI-bearing fields, 310 unique DOIs (unchanged; no new
  citations were introduced).
  holdout_guard PASS — 3/3 invariants, re-run after the refit and after regeneration.
  fit_target_gate PASS — both checks. It now enumerates
  `matrix_observability_refit_pratap_singh.json: leverage=per_row_recovery, 2 target(s)`,
  and the per-row fit-target set is
  `['cys_ribose_140C_Hofmann1998', 'pea_isolate_40C_PratapSingh2021', 'soy_isolate_40C_PratapSingh2021']`.

  FINAL SUITE: see the line appended below.

### [P] CARRIED FORWARD

  1. The 1-hexanol observability factors are composed entirely of fabricated values
     (`0.143 / 0.063` = 120 ppb / 80 ppb, neither in the paper) and have no anchor to refit
     against. Retire, or find a real soy/pea 1-hexanol pair.
  2. The compound-specific matrix calibration registry is unreachable on the
     `matrix_precursor_augmented` path (SMILES-vs-name lookup). Either label the species or
     canonicalise the lookup — until then the internal snapshots cannot detect drift in it.
  3. The pea/soy CLASS-level aldehyde anchors were not moved with the refit. They remain
     correct as a RATIO but their absolute scale is now inconsistent with the
     compound-specific hexanal lane.
  4. Bi 2020 (1260 ppb) vs Liu 2023 (15-180 ppb band) disagree by 24x on nominally the same
     system and no observability factor satisfies both. Resolving which is representative of
     commercial PPI is the single highest-value retrieval for this lane, ahead of Trikusuma.
  5. Trikusuma 2020 remains the last content-unverified pillar of the matrix lane and the
     only fit-recovery benchmark still scoring a full `pass`.

### FINAL SUITE (Wave O)

  **1245 passed, 1 skipped, 2 xfailed, 0 FAILED** in 833.90 s
  (`tests/unit tests/scientific tests/integration tests/scripts`, documented conda path).
  Entry state was 1242/0 (Wave M) + the 3 new tests in
  `tests/scientific/test_matrix_observability_refit_wave_o.py` = 1245. The 1 skip and both
  xfails are the declared, strict-marked ones from Wave J2. All three gates re-run green
  after the suite. NOT COMMITTED, NOT STASHED — handed to the orchestrator as instructed.

---

## Wave Q — verification of `data/lit/timeseries/` (2026-08-27)

Scope: verify the four harvested kinetic time-series files against their source papers,
point by point, so they can be committed and later used for calibration. They were produced
by an agent killed before verifying anything; every number was treated as unverified.
Files touched: `data/lit/timeseries/*` only. Nothing committed. Nothing wired into the model.

### Verdict table

| File | Series | Points | Status | Evidence |
|---|---:|---:|---|---|
| `brands_sugar_casein_120C_pH68.yml` | 20 | 60 | **verified_table** | Brands (2002) WUR thesis retrieved (edepot.wur.nl/199005). Tables 6.1 (PDF p.105), 6.2 (p.107), 3.1 (p.44) rendered at 4x and re-read as images — the OCR text layer of a bitonal scan is not trustworthy on its own. Every digit, every ±, every SEM matched. |
| `martins2003_DFG_amadori_degradation.yml` | 8 + table | 106 | **verified_table + verified_figure_shape** | Martins (2003) WUR thesis retrieved (edepot.wur.nl/121418). Table 4.1.1: 32 values + 12 `n.a.` cells exact. Figs 4.1.1/4.1.2: all 76 markers re-extracted independently, reproduced to ≤0.01 mmol/L. |
| `martins2005_glucose_glycine_100C_pH68.yml` | 10 | 183 | **verified_figure_shape** | Thesis Fig 5.9 (PDF p.126) re-extracted from the raw content stream. **All 183 points reproduced exactly**, to the last printed digit. Marker census 24/species = 8 times × 3 replicates. |
| `martins2005_glucose_glycine_80_100_120C_pH68.yml` | 42 | 418 | **corrected** | Thesis Fig 5.10 (PDF p.131). 415/418 reproduced; **3 phantom points removed** (below). Nine of ten species matched perfectly. |

All four sources were **actually retrieved**. The thesis route paid off: both Wageningen
theses download openly and both are the primary form of the published chapters.

### Citations

All six DOIs resolve on CrossRef with exact title/journal/volume/page/author agreement:
`10.1016/j.foodchem.2004.04.006`, `10.1016/s0008-6215(03)00173-3`,
`10.1016/s0008-6215(03)00174-5`, `10.1021/jf9907586`, `10.1021/jf010789c`,
`10.1021/jf001430b`. **Zero citation contamination in this directory** — a contrast with the
30–45% rate the 2026-08 audit found elsewhere in the reference set.

**The Brands file does not need a DOI found.** Its source is a PhD thesis, which has none;
it is openly downloadable. The two chapters supplying every number are published articles
whose DOIs were already present and now verified (JAFC 2000 mutagenicity → Tables 6.1/6.2;
JAFC 2002 melanoidin quantification → Table 3.1). The task brief's guess that these data come
from Brands & van Boekel 2001 JAFC `10.1021/jf001430b` is **wrong**: that is a different
chapter of the same thesis and none of its data are in the file. It is the right *target*
for the next retrieval, not the source of anything already here.

### The one data defect found (in 767 values)

`martins2005_glucose_glycine_80_100_120C_pH68.yml`, glycine 120 °C — three values that do
not exist in the figure, now removed:

- `[10, 196.1]` — at t=10 the figure draws exactly 3 markers: two coincident at 190.6, one
  at 191.9. No third position. 196.1 is the t=5 value.
- `[15, 179.5]` — at t=15 the figure draws 9 markers (3 T × 3 replicates) at 6 distinct
  positions; none is at 179.5. 179.5 is the t=20 value.
- `[20, 185.1]` — at t=20 only 120 °C is sampled; 3 markers, two at 178.2, one at 179.5.
  185.1 is the t=15 value.

Each removed value duplicates an adjacent time point — the series had been padded so every
time point showed three replicates. Detectable only because the re-extraction recorded
marker **multiplicity** per position: a position with zero drawn markers cannot be data.
`n_points` 15 → 12. Source of correction: content-stream enumeration of Fig 5.10 panel B,
thesis PDF p.131.

### Other corrections

Data added: 6 values (Brands Table 3.1 `a/b'` microanalysis column, previously skipped).

Metadata: Brands ISBN `90-5808-579-4` → `90-5808-591-0`; unsupported "137 pp" removed
(133-page scan ending at printed p.127, no page count stated anywhere); Brands Table 6.1
p.101 → p.99 and Table 3.1 p.36 → p.37; scan "200 dpi" → measured ≈160 dpi 1-bit; Martins
ISBN `90-5808-923-4` → `90-5808-823-5` (both files); DFG file's "filled-diamond series"
exclusion note rewritten (no diamonds exist on that page — the fits are stroked polylines,
the exclusion was right but its description was not); Fig 5.10 dotted fit curve relocated
from panel B to panel C; Fig 5.10 read-off agreement restated from "under 1% of full scale"
to the measured worst cases (2.2% for DFG at 60–120 min).

### Verified and vindicated

- The DFG file's flag that the thesis M&M is internally inconsistent (`0.237 g, 10 mmol`;
  0.237 g of DFG = 1.0 mmol) is **confirmed verbatim**. The error is the thesis's, and it was
  correctly flagged rather than silently fixed.
- Every `system:` block (concentrations, buffer, pH, vessel, volume, filtration, replication)
  confirmed verbatim against the retrieved Materials and Methods of both theses.
- Fig 5.9's formic-vs-acetic assignment — which the caption cannot settle — confirmed from
  the thesis body: *"acetic acid was always formed in higher concentrations than formic
  acid"* (printed p.109) and *"the low amount of formic acid detected relatively to acetic
  acid at pH 6.8"* (p.122).
- Melanoidin ε = 0.64 ± 0.03 L·mmol⁻¹·cm⁻¹ at 470 nm confirmed (printed pp.37, 52).
- The Fig 5.10 temperature attribution — the file's most reconstructed claim — holds. The
  100 °C series agrees with the independent Fig 5.9 extraction at every shared time point for
  all ten species; 5/10/20 min are 120 °C-only and 150/180 min are 80 °C-only, so those need
  no splitting; the genuinely ambiguous points were already parked under
  `ambiguous_80_or_100` and that parking is correct.
- Marker census: exactly 66 markers per species in Fig 5.10 (= 22 time-temperature cells × 3
  replicates) and 24 per species in Fig 5.9 (= 8 times × 3), so replication was 3-fold, not
  merely the "at least duplicate" the M&M claims.

### What remains unverifiable, and what a human should pull

Nothing in these files is now unverifiable. What is **absent** and worth retrieving, in order:

1. **Brands & van Boekel 2001 JAFC `10.1021/jf001430b`** — the glucose/fructose-casein
   multiresponse trajectories. In the thesis they are figures in a 1-bit ≈160 dpi scan with
   zero vector content, so they are genuinely unreadable there; reading points off them would
   be eyeball estimation and was correctly not attempted. Someone with ACS access should
   check whether the publisher PDF is vector. **This is the biggest gap** — the Brands file
   currently holds endpoints (browning, Ames, C/N), not trajectories.
2. **Martins thesis Figs 4.1.3–4.1.7** (PDF pp.78–82) and the Ch. 4.2 fit figures
   (pp.93–105): deoxyosones, MGO, sugars, acids, melanoidins during DFG degradation at four
   conditions. Vector, same PDF, same method. Skipped for time only.
3. **90 °C and 110 °C glucose/glycine series** — run and fitted by the authors, plotted
   nowhere. Needs raw data from Wageningen. Would give a five-point Arrhenius set.
4. **Underlying tables for Figs 5.9/5.10** — never published in either the thesis or the
   Food Chemistry paper.

### Standing caveat

Three of four files are **figure read-offs**. The re-extraction proves the numbers faithfully
recover where the authors' plotting program placed each marker; it cannot prove those
positions equal the authors' measurements. The ±0.3% of panel full scale is a **floor**, not
a total. For Fig 5.10 there is a second layer: the temperature attribution is a
reconstruction, corroborated but not labelled by the figure. Both caveats are recorded in
each file's `residual_uncertainty` block and in `data/lit/timeseries/README.md`.

`data/lit/timeseries/README.md` written: directory purpose, per-file status table, the
fabrication warning pointing at `AUDIT.md`, the explicit statement that nothing here is wired
into calibration, the verification method, and the retrieval backlog above.

## Bi-vs-Liu adjudication (2026-08-27, orchestrator, owner-requested)

Wave O's sharpest [P] — "Bi 2020 and Liu 2023 disagree by 24x on the pea hexanal lane
and no constant can satisfy both" — is RESOLVED, and the resolution is that the repo's
Liu numbers are untraceable to their source. Primary evidence retrieved:

1. The "Liu, Y. (2023 thesis)" anchor is Yaozheng Liu, "Flavor Chemistry of Pea
   Proteins", NC State MS thesis (2021; the derived paper is Liu, Cadwallader & Drake,
   Food Chemistry 406:134998, 2023, DOI verified). FULL THESIS RETRIEVED (216 pp, via
   repository.lib.ncsu.edu). Method: 24 commercial pea proteins (9 concentrates, 15
   isolates), rehydrated to 10% solids; quantification by HS-SPME-GC-MS/MS external
   standard curves with response factors (R^2>0.95), units ug/L of the 10% slurry
   (Table 2.7). CAVEAT: curves built in DI WATER, not protein matrix, so protein
   binding is uncorrected -> values are, if anything, LOWER bounds on total hexanal.
2. Table 2.7 hexanal across the 9 quantified samples: 2445 to 52454 ug/L (typical
   2400-12000). Nonanal: 0.188 to 3.42 ug/L. 2-pentylfuran: 25.5 to 1628 ug/L.
3. THE REPO'S LIU BAND (hexanal 15-180 ppb; nonanal band mid 15.81) MATCHES NOTHING
   in the thesis: hexanal is ~50-300x HIGHER than our band; nonanal ~5-80x LOWER.
   The 15-180 band appears confabulated or transcribed from an unidentified source.
4. Pratap-Singh 2021 (verified earlier): standards spiked INTO the 1:7 w/v protein
   slurry (matrix-matched, hexanal-d12 IS) -> per-slurry basis, 1138/1622 ug/L.
   COMPARABLE BASIS to Liu (10% vs 12.5% solids).

VERDICT: on comparable slurry bases the two believable measurements AGREE within
~2-5x at the low end of Liu's inter-lot range: Pratap-Singh pea PPI 1138 ug/L; Liu
commercial pea proteins 2445-52454 ug/L. Commercial pea-protein slurries carry
hexanal in the THOUSANDS of ug/L. Bi 2020 (raw/roasted PEA FLOUR, not isolate;
1260 = OAV 280 x our own 4.5 ppb ODT, i.e. derived not measured) is consistent at
the low end for a different matrix. The "24x literature conflict" was a conflict
between a real measurement and a phantom band. TRUST ORDER for the pea lane:
Pratap-Singh (matrix-matched SIDA-style) >= Liu thesis Table 2.7 (water-calibrated,
lower bound, huge real inter-lot spread) >> Bi (derived, wrong matrix) >>> the
retired 15-180 band (no identifiable source).

QUEUED (after Wave P lands - its regeneration reads these files): correct
external_validation_liu_2023_ppi_offnote_baseline.{json,yaml} + payload anchors to
Table 2.7 bands (hexanal band 2445-52454, geometric mid ~11330; nonanal band
0.188-3.42, geometric mid ~0.80; provenance thesis_table, water-calibration caveat),
re-derive ranks, re-score hold-out, sync docs/tests. Note Wave O's post-refit
prediction 1013 ug/L for liu-hexanal lands 1.2x UNDER Liu's lowest lot vs 11x under
the geometric mid - the refit direction is VINDICATED by the corrected target.
[P] the nonanal correction will worsen that point (model 171.7 vs ~0.80 mid) -
consistent with the oleate-nonanal topology fix Wave P is landing (nonanal should
mostly vanish when routed off linoleate).

## Wave P — evidence-grounded chemistry additions (2026-08-27, owner-approved)

Six owner-approved items. Every step atom-balanced and RDKit-verified, every constant either
literature-anchored or labelled ESTIMATED/UNCONSTRAINED, every DOI CrossRef-verified BEFORE it
was written anywhere. Nothing was tuned to reproduce a prior output, and the two places where
the numbers got worse are reported as such.

ORDER OF EXECUTION, and it matters. The items were EXECUTED 4, 6, 2, 3, 5, then 1 — not in the
order they were listed. Item 1 is a refit, and Wave L1's own finding about Wave H was that
fitting a barrier before a route change fits the wrong route more precisely. The same logic
applies inside a wave: items 2 and 6 both move the ribose/cysteine system, so the refit was run
LAST, against the network that actually ships. Item 7 needed no code change (see below).

### (a) THE FINDING OF THIS WAVE — parallel channels cannot add

Item 2 added the MFT channel Hofmann & Schieberle 1998 measured as their HIGHEST-yielding
system. Measured contribution to the shipped prediction: **EXACTLY ZERO**.

  On cys_ribose_140C_Hofmann1998, at the shipped constants, 2026-08-27:
      pentodiulose lane alone   242.38 ppb
      C2+C3 lane alone           71.02 ppb
      BOTH lanes                242.38 ppb      <- the max, not the sum

`src/recommend.py::predict_from_steps` relaxes to the LOWEST-SPAN path per product (see the
`if path_span < current_span` branch); it does not sum flux over parallel channels. Had the two
added, MFT would read 313.39 ppb = 1.09x under instead of 1.41x under. Consequences, stated
plainly:
  * every compound-level claim in this repository rests on ONE route, whichever is fastest;
  * enriching the network CANNOT improve a number unless the new route is faster than the old;
  * "the model has N routes to X" is not a statement about the prediction of X.
REPORTED, NOT FIXED. Making the propagator additive is a model-wide recalibration event with
its own regeneration and belongs to its own wave. Pinned in
`tests/scientific/test_wave_p_chemistry_2026_08.py::test_the_second_mft_channel_contributes_exactly_zero_to_the_prediction`,
which fails the day someone changes it — and that failure is the notification.

### (b) A NAMING HAZARD, found while placing the new families

`src/conditions.py` selects the pH-ionisation and Labuza water-activity corrections by
SUBSTRING MATCH on the family name ("strecker", "amadori", "schiff", "pyrazine", "furfural").
The natural name for item 6's new family was `Furanone_Strecker_Reduction` — which would have
silently collected BOTH corrections and suppressed DMHF by ~480x at pH 5.5 / aw 0.95, as a side
effect of a spelling. Its sibling `Furanone_Formation`, the same chemistry, receives neither.
The family is therefore called `Furanone_Amino_Acid_Reduction`, and a test now pins that none
of the six Wave P families picks up either correction. THE UNDERLYING COUPLING IS UNTOUCHED and
is an open item: it means a family's name is a load-bearing calibration input. [P]

### (c) ITEM 1 — refit `thiol_addition_pentodiulose` against Hofmann1998

  Generator: `scripts/generators/refit_thiol_addition_pentodiulose_hofmann.py` (mirrors
  `refit_sulfur_barriers_hofmann.py`'s honest-profile style: same objective, same
  identifiability / materiality / conservative-edge rules, same grid step, plus an explicit
  BOUNDARY-HONESTY check the older script lacked).
  Record: `results/validation/sulfur_barrier_refit_pentodiulose.{json,md}` (un-ignored in
  .gitignore). It SUPERSEDES `sulfur_barrier_refit_hofmann.{json,md}`, which profiles a family
  no step emits; the old record is left in place for provenance and still declares the same fit
  target, so the fit-target accounting is unchanged by the new file's arrival.

  ADOPTED: 28.60 -> **26.35** kcal/mol.
  * range [23.30, 29.65] = thiohemiacetal_formation .. thiol_addition_hexose, both already in
    FAST_BARRIERS. Same range the retired knob used, deliberately, so the two fits compare.
  * ONE free parameter. The Wave H script fitted FOUR barriers against these two rows (2.0
    parameters/row); this is 0.5, exactly at FIT_LEVERAGE_THRESHOLD, so the benchmark stays
    `per_row_recovery` and out of the honest coverage numerator AND denominator either way.
    What changed is how much freedom the fit had, not the benchmark's status.
  * objective (mean |log10| over the 2 rows) **0.2192 -> 0.0935 dex**; profile span 0.3224 dex.
  * **THE PROFILE MINIMUM IS AT THE RANGE FLOOR** (argmin 23.30, `hit_a_bound = true`). The
    data want a value the table's own sulfur-addition class envelope does not contain: the
    profile SATURATES, and with the knob pinned at the floor the objective is still 0.0836 dex.
    26.35 is the conservative edge of the indifference band, NOT the optimum. The residual is
    not removable by this barrier.
  * MFT 151.87 -> **242.38** ppb vs 342: 2.2519x under -> **1.4110x under**.
  * FFT 243.72 -> **217.99** ppb vs 200: 1.2186x over -> **1.0900x over**. **FFT CO-MOVES, AND
    IT WAS NOT FITTED.** Over the search range MFT spans 3.53x and FFT spans 1.25x, in
    OPPOSITE directions: lowering the MFT barrier diverts shared upstream sugar flux into the
    MFT lane and FFT falls. The two rows are therefore not independent evidence — which is
    precisely why ONE knob was fitted against them and not two.
  * REPORTED, NOT FITTED (objective span over +/-2 kcal, at the adopted point):
        furanone_cyclisation      0.0000 dex  (exactly zero derivative here)
        deoxyosone_reduction      0.3850 dex  (in SERIES with the fitted step; fitting both
                                               would split one measurable quantity in two)
        thiohemiacetal_formation  0.0033 dex  (FFT-side)
        thiol_dehydration         0.0285 dex  (FFT-side)
        thiol_addition            0.1011 dex  (labels only the demoted lumps)
  * THE CAVEAT, CARRIED VERBATIM. The fit target's own `content_verification_note` (Wave K)
    says the 342/200 ppb require a mol%->ppb conversion that is "NOT documented anywhere in
    this repo" and could not be confirmed against the paper. That sentence is now quoted
    VERBATIM inside the constant's FAST_BARRIERS rationale, inside the JSON record, and in the
    generator, with the consequence spelled out: the constant is fitted against an anchor whose
    unit conversion is UNVERIFIED, and if the conversion is wrong THE ERROR IS LOCALISED IN THIS
    ONE BARRIER rather than distributed across the route. The generator ASSERTS the note is
    still present and unchanged before it runs, so the caveat cannot go stale silently.
  * The benchmark's own contract is UNTOUCHED and still FAILS — max_ratio 1.4110 is now INSIDE
    its 1.45 threshold but MALE 0.0935 is outside 0.09, i.e. it fails on one criterion instead
    of two.
  * Note 26.35 vs the retired 26.85 that Wave H fitted THROUGH the contradicted route. The two
    fits are independent (different route, same benchmark, same rules) and the near-agreement is
    a coincidence of the shared upstream trunk, NOT evidence for the retired route. A test
    asserts the adopted value is not 26.85.

### (d) ITEM 2 — the C2 + C3 recombination lane to MFT

  Hofmann & Schieberle 1998 (10.1021/jf9705983, abstract retrieved VERBATIM via Europe PMC,
  PMID 10554225 — CrossRef-verified): "the highest yields for MFT (1.4 mol %) when
  hydroxyacetaldehyde and mercapto-2-propanone were reacted for 6 min at 180 degrees C in the
  absence of water". The network could not express this at all: mercapto-2-propanone was not a
  species. Wave L1 listed it as one of the top three completeness gaps.

  THREE STEPS, all exact and RDKit-verified:
      pyruvaldehyde + H2S + 2[H] -> 1-mercapto-2-propanone + H2O
          C3H4O2 + H2S + H2 -> C3H6OS + H2O          both sides C3H8O2S
          family `Mercaptoketone_Formation`
      glycolaldehyde + 1-mercapto-2-propanone -> 4,5-dihydroxy-3-mercapto-2-pentanone
          C2H4O2 + C3H6OS -> C5H10O3S                both sides C5H10O3S
          family `Mercaptoketone_Aldol_Addition`  (ADDITION, not condensation — the named
          literature intermediate retains both oxygens)
      4,5-dihydroxy-3-mercapto-2-pentanone -> MFT + 2 H2O
          C5H10O3S -> C5H6OS + 2 H2O                 both sides C5H10O3S
          family `Mercaptoketone_Cyclodehydration`
  Overall C2H4O2 + C3H6OS -> C5H6OS + 2 H2O, exact, with NO reducing token on the
  recombination itself; the token is consumed once, upstream, in step 1.

  MECHANISM CITATION. The 1998 abstract states no mechanism (it is a yield paper; the body is
  paywalled and was NOT retrieved). Cerny 2015, "The role of sulfur chemistry in thermal
  generation of aroma", 10.1016/b978-1-78242-103-0.00009-6 (CrossRef-verified; FULL TEXT
  retrieved) renders their Figure 9.6 verbatim: "The formation of 2-methyl-3-furanthiol (16)
  starts with aldol reaction of hydroxyacetaldehyde (18) and mercaptopropanone (19) to give
  4,5-dihydroxy-3-mercapto-2-pentanone (20). The intermediate 20 cyclises and dehydrates to
  yield 2-methyl-3-furanthiol (16)." Both intermediates are NAMED there; neither is invented.

  WHERE THE C3 PARTNER COMES FROM — AND THE ROUTE THAT WAS REJECTED. Two candidates were
  balanced and researched:
      (a) hydroxyacetone + H2S -> 1-mercapto-2-propanone + H2O   [EXACT, no token]
      (b) pyruvaldehyde + H2S + 2[H] -> 1-mercapto-2-propanone + H2O
  (a) is the CHEAPER step and needs no reducing equivalent. It was REJECTED: a literature search
  found NO support for direct thiol-for-hydroxyl exchange on a simple alpha-hydroxyketone under
  Maillard conditions. (b) is what the literature states, three separate times in Cerny 2015:
  "An alternative formation pathway to alpha-mercaptoketones passes through dicarbonyl
  compounds ... These can add hydrogen sulfide ... and then undergo reduction. Actually,
  3-mercapto-2-butanone was found in the model reaction of hydrogen sulfide with
  2,3-butanedione, and 3-mercapto-2-pentanone and 2-mercapto-3-pentanone in the reaction with
  2,3-pentanedione (Mottram et al., 1995)". Primary model-system evidence: Mottram, Madruga &
  Whitfield 1995, 10.1021/jf00049a035 (CrossRef-verified, JAFC 43:189-193).
  HONEST LIMIT: the demonstrated cases are the C4 and C5 dicarbonyls. The C3 member is the
  analogous instance of a mechanism stated GENERALLY for "dicarbonyl compounds" and demonstrated
  for its two nearest homologues. That is written into the docstring and the barrier rationale.
  1-mercapto-2-propanone is itself reported as a cysteine/ribose and cysteine/glucose Maillard
  product in Cerny 2015's own product lists (attributed there to Hofmann & Schieberle 1995) —
  but NOT in the 1995 abstract, which was retrieved and does not name it. Flagged.

  NOT MOISTURE-GATED, and why. The 1.4 mol % was measured "in the absence of water" and
  `ReactionConditions` does carry `water_activity`, so a gate is mechanically possible. It was
  declined: a single dry-condition yield does not establish the channel is ABSENT in water, and
  Cerny 2015 says both intermediates "are known to occur in heated L-cysteine/glucose systems",
  which are aqueous. Worse, the existing aw machinery is only reachable through the family-name
  substring matcher, i.e. by choosing a name containing "furfural" — a naming trick, not a
  model. The channel competes on barriers alone and THE MODEL DOES NOT EXPRESS THE MOISTURE
  DEPENDENCE THE MEASUREMENT WAS MADE UNDER. Stated as a limitation. [P]

  REACHABILITY, MEASURED — AND THE MEASUREMENT CORRECTED A GUESS I HAD ALREADY WRITTEN DOWN.
  The lane needs glycolaldehyde, which `_retro_aldol_fragmentation` emits only from the PENTOSE
  deoxyosone, so the docstring first claimed the lane was pentose-only. Enumeration says
  otherwise: in glucose+cysteine the identical species arrives via glyoxal -> Strecker ->
  2-aminoethanal -> hydrolytic deamination, because pool identity is by canonical SMILES and not
  by label. The docstring was corrected to the measurement. While there,
  `_deamination_step._known_ketol_labels` gained `"O=CCO": "glycolaldehyde"` so the molecule
  stops appearing as `deaminated-2-aminoethanal-ketol` — a READABILITY fix with zero numeric
  effect (verified: the panel is byte-identical across it).

### (e) ITEM 3 — norfuraneol -> 2,3-pentanedione -> 2-mercapto-3-pentanone

  Wave N retired norfuraneol from the MFT lane and left it with ZERO consumers. Cerny & Davidek
  2003 (10.1021/jf026123f, abstract VERBATIM, PMID 12696962) name what it actually does:
  "Whereas 2-mercapto-3-pentanone was found unlabeled and hence originated from
  4-hydroxy-5-methyl-3(2H)-furanone, its isomer 3-mercapto-2-pentanone was formed from both
  4-hydroxy-5-methyl-3(2H)-furanone and ribose."

  THE ISOMER MATTERS AND THE LANE EMITS ONLY ONE. 2-mercapto-3-pentanone is assigned
  exclusively to norfuraneol; 3-mercapto-2-pentanone has a second, ribose-derived origin (the
  1,4-dideoxyosone route already in the network), so emitting it here would double-count that
  flux. A test asserts it is NOT emitted.

  THE MECHANISM IS NOT IN THAT PAPER — it proposes a pathway only for MFT and
  3-mercapto-2-pentanone. Whitfield & Mottram 1999 (10.1021/jf980980v, CrossRef-verified, JAFC
  47:1626-1634, abstract VERBATIM, PMID 10564029) reacted norfuraneol with cysteine or H2S at
  pH 4.5 / 140 C / 60 min: "The main non-sulfur compounds were 2,3-pentanedione,
  2,4-pentanedione, and 3,4-hexanedione." So the ring opens to the C5 alpha-dicarbonyl FIRST and
  sulfur enters afterwards, by the same alpha-mercaptoketone mechanism as item 2.

  TWO STEPS, exact and RDKit-verified:
      norfuraneol + 2 x 2[H] -> 2,3-pentanedione + H2O
          C5H6O3 + 2 H2 -> C5H8O2 + H2O              both sides C5H10O3
          family `Furanone_Reductive_Opening`
      2,3-pentanedione + H2S + 2[H] -> 2-mercapto-3-pentanone + H2O
          C5H8O2 + H2S + H2 -> C5H10OS + H2O         both sides C5H12O2S
          family `Mercaptoketone_Formation`  (SAME family as item 2 — same mechanism, and
          Mottram 1995 demonstrated it on this exact substrate)

  THE ONE-STEP LUMP WAS WRITTEN OUT AND REJECTED. `norfuraneol + H2S -> 2-mercapto-3-pentanone
  + 2 H2O` needs SIX hydrogens (C5H8O3S vs C5H14O3S), i.e. THREE reducing-equivalent tokens on
  a single step. NOTE: the brief for this wave stated the requirement as +4H / two tokens; the
  arithmetic is +6H / three, verified independently twice. Three tokens on one step is not an
  elementary step by any reading; the two-step literature route accounts for the same three
  across two steps with a MEASURED intermediate in between.

  SENSORY LAYER: 2-mercapto-3-pentanone is quantified in heated meat by Kerscher & Grosch 1998
  (10.1021/jf970892v, CrossRef-verified: beef 20-44 ug/kg, pork 11-14, lamb 10, chicken 13) but
  NO PEER-REVIEWED ODOUR THRESHOLD FOR IT COULD BE RETRIEVED — it is absent from Cerny 2015's
  threshold tables and from Blank's sulfur-odorant table, both of which carry the ISOMER
  (3-mercapto-2-pentanone, "sulfury, catty", 0.7 ug/kg in water). It is therefore NOT registered
  as a scored sensory target: adding a target with a fabricated threshold is the exact failure
  mode this campaign exists to remove. It ships as a network species and a sulfur sink only. [P]

### (f) ITEM 4 — nonanal comes off the OLEATE pool

  `LipidProfile.oleic_acid_pct` was declared, populated for pea and soy, carried through the
  calibration registry, and READ BY NO CODE PATH IN THE REPO. Nonanal is in NEITHER linoleate
  hydroperoxide isomer's product list (Miyazaki 2023, 10.1093/bbb/zbac189, read in FULL TEXT):
  13-HpODE -> hexanal, vinyl hexanoate, 2-pentylfuran, 4,5-epoxy-2-decenal; 9-HpODE ->
  2,4-decadienal, octanoic acid, 2-octenal, hexanal. Nonanal is the C9 fragment of the OLEATE
  double bond.

  CITATION CORRECTION MADE WHILE HERE. `docs/validation/isotope_topology_evidence.md` attributes
  10.1021/jp0500900 to "Hearn & Smith (2007)". CrossRef resolves that DOI to **Hung, Katrib &
  Martin (2005)**, J. Phys. Chem. A — the dossier has the TITLE right and the authors/year
  wrong. The abstract (Europe PMC, PMID 16833788) does carry the quoted figure verbatim:
  "1-nonanal (30 +/- 3% carbon yield)". The correct attribution is used in the code; the dossier
  is a Wave L1 deliverable and was NOT edited (it is on this wave's do-not-touch list). The
  companion Reynolds citation is 2006, not 2007 (the dossier says both in different places). [P]

  IMPLEMENTED: `predict_hexanal_generation` now returns `total_hydroperoxide_oleate` alongside
  `total_hydroperoxide`, computed with the SAME kinetics from `oleic_acid_pct`; nonanal is
  scaled off it. `MARKER_HYDROPEROXIDE_POOL` / `hydroperoxide_pool_key_for_marker` are the
  single source of truth and are consumed by both matrix paths in `benchmark_validation.py` and
  by `derive_matrix_sigma_from_residuals.py`. The BRANCHING RATIO (0.15) was NOT touched — only
  the pool it multiplies.

  KNOWN, DELIBERATELY UNIMPLEMENTED: the oleate pool uses the same rate constant as the
  linoleate pool, which is generous to nonanal (monoenoic C18:1 oxidises ~an order of magnitude
  slower than dienoic C18:2). Applying that ratio would push predicted nonanal down another
  ~10x, i.e. FURTHER TOWARDS the hold-out — which is exactly why it was not done in the same
  pass. Predicted nonanal is expected to remain biased HIGH. [P]

  OLD -> NEW, EVERY ROW THAT MOVED:
    EXTERNAL HOLD-OUT (n=200, seed 0; hold-out never fitted, holdout_guard re-run PASS)
      | point                          | measured | p50 old  | p50 new | fold old  | fold new |
      | li_2026_hme / nonanal          |    72.66 | 19809.0  |  8596.4 | 272.6257x | 118.3093x|
      | liu_2023_ppi / nonanal         |    15.81 |   171.70 |   75.55 |  10.8603x |   4.7785x|
      the other SIX points are byte-identical.
      Each moved by exactly its matrix's oleic/linoleic ratio (soy 0.4340, pea 0.4400).
      AND THE HEADLINE DID NOT MOVE AT ALL: median |log10| 1.6296 dex / median fold 42.6159x
      UNCHANGED, ci_coverage_hits 4/8 UNCHANGED, max_fold_error 2474.3839x UNCHANGED,
      pre-widening genuine_extrapolation_hits 1/5 UNCHANGED. The median sits between two points
      that did not move, and the two that did were outside the interval before and after.
      A real mechanistic correction improved a quarter of the hold-out by 2.3x each and moved
      the published number by nothing. That is the clearest argument in the repo for reading the
      per-point table alongside the headline.
    IN-PANEL
      pea_isolate_uht_140C_Trikusuma2019 nonanal 24.00 -> 10.56 ppb = 2.2727x under, which is
      EXACTLY 1/(22.0/50.0). Benchmark `pass` -> `scale-gap`. THE OBSERVABILITY FACTOR WAS NOT
      REFITTED: doing so would re-absorb the correction into the same constant and make the fix
      invisible, and Trikusuma 2019 is still the last content-unverified pillar of the matrix
      lane, so there is no verified anchor to refit against. Dated note added to
      `src/matrix_calibration_registry.py`.
      -> THIS BREAKS THE LAST FIT-RECOVERY PASS IN THE PANEL: fit-recovery passes 1/4 -> 0/4,
      aggregate (strict `pass`) 5/14 -> 4/14. Every remaining pass is an internal synthetic row.
      The four Internal2026/ProtocolPilot2026 snapshots' Nonanal rows moved by the same ratio.
    MATRIX SIGMA (`derive_matrix_sigma_from_residuals.py`, leave-lane-out)
      One of the five scored residuals is the Trikusuma nonanal row: uncalibrated prediction
      3238.93 -> 1425.13 ppb. rms_ln_sigma 3.2520 -> **3.0166**; bias_fold 3.9065 -> 3.3149;
      centered_sd 3.3013 -> 3.0951; 90% CI [2.1856, 6.7958] -> **[2.0274, 6.3038]**.
      The shipped 2.86 was NOT moved and is still inside. It moved TOWARDS the shipped value;
      that is a consequence of correcting a substrate, not a reason to narrow the prior.

### (g) ITEM 5 — fructose reaches HMF by its own route

  The model sent fructose Heyns -> hexose 3-deoxyosone -> HMF, i.e. through the precursor an
  isotope study excludes for exactly this sugar. Perez Locas & Yaylayan 2008 (10.1021/jf8010245,
  CrossRef-verified, abstract VERBATIM, PMID 18611024): "both fructose and sucrose showed much
  higher conversion rates than 3-deoxyglucosone thus precluding it as a major precursor of HMF
  in fructose and sucrose solutions ... sucrose degrades into glucose and a very reactive
  fructofuranosyl cation."
  CAREFUL WITH THAT CITATION, and the care is recorded in the code: its title and design are
  about SUCROSE and the cation it proposes arises from GLYCOSIDIC CLEAVAGE. It is cited for
  exactly two things — the 3-DG exclusion for fructose, and the cation's identity. The
  ring-retained route from FREE fructose is Antal, Mok & Richards 1990
  (10.1016/0008-6215(90)84096-d, CrossRef-verified, Carbohydr. Res. 199:91-109, abstract
  VERBATIM), which reviews both hypotheses and finds for the ring-retained one.

  IMPLEMENTED, minimally: for a KETOSE the 1,2-enolisation HMF step is replaced by
      D-fructose C6H12O6 -> HMF C6H6O3 + 3 H2O      [EXACT, family `Fructofuranosyl_Dehydration`]
  The Heyns product and the 3-deoxyosone are still formed (they feed the retro-aldol and
  Strecker lanes); only the claim that fructose's HMF comes through them is removed. In a system
  with both sugars the glucose limb still emits deoxyosone -> HMF and pool de-duplication means
  it appears once. NO INTERMEDIATE WAS INVENTED: Amarasekara 2008 (10.1016/j.carres.2008.09.008)
  names a structurally defined one, but in DMSO with the solvent as catalyst, so it is recorded
  and not used. NO RATE ADVANTAGE ENCODED: the barrier is set EQUAL to the `dehydration` class
  value the replaced step carried (28.0), even though the paper measured fructose converting
  FASTER. Conservative direction, deliberately not tuned.

  MEASURED IMPACT ON SCORED ROWS: **ZERO**. The only benchmarks containing fructose are the two
  QUARANTINED PMC9905368 files. This is a topology repair with no scoreboard consequence, and it
  is reported as one rather than left unmentioned because it moved nothing.

### (h) ITEM 6 — the furaneol `[HH]` pool gate is gone

  Red-team H4's second half, and Wave L1's topology risk 3. The hexose 1-deoxyosone -> DMHF step
  was pool-gated on the `[HH]` token whose only producer in a cysteine-free system is the
  pyrazine aromatisation, so predicted furaneol from glucose was a downstream dependent of
  pyrazine chemistry — while Wang & Ho 2008 (10.1021/jf8012025) CONFIRM the route with CAMOLA.

  A `_sugar_reductant_pool` WAS INVESTIGATED AND DECLINED. Research verdict, and it is the
  honest answer to the question the brief asked: the literature establishes that Maillard
  reductones ACT as reducing agents — including the 1-deoxyosone itself, measured (Kanzler,
  Haase & Kroh 2014, 10.1021/jf404322r; Kanzler et al. 2016, 10.1021/acs.jafc.6b03398) — but it
  does NOT pin a specific balanced couple with both the reductone and its dehydroreductone named
  and structurally defined: Kanzler 2016 finds the oxidised partner is a TRANSIENT TRICARBONYL
  that fragments to lactic/glycolic/glyceric acid, so what was characterised is the fragments,
  not the couple. Writing that couple would have meant inventing a species purely to carry
  bookkeeping — the exact move Wave I declined for the same reason. A search for a
  "reductone disproportionation" claim returned zero primary retrievals and it is NOT asserted.

  WHAT WAS DONE INSTEAD IS BETTER: the token was not re-sourced, it was REMOVED, because the
  accepted mechanism NAMES the reductant and it is the amino acid. Blank & Fay 1996
  (10.1021/jf950439o, CrossRef-verified, abstract VERBATIM): the mechanism is "based on
  decomposition of the Amadori compound via 2,3-enolization, chain elongation by the Strecker
  aldehydes, and REDUCTION of the resulting acetylformoin-type intermediates to the target
  molecules". Kerler, Winkel, Davidek & Blank 2010 (10.1002/9781444317770.ch3, CrossRef-verified,
  FULL TEXT) name the donors: "acetylformoin ... is an effective precursor for
  4-hydroxy-2,5-dimethyl-3(2H)-furanone (Furaneol). The amounts of Furaneol obtained from
  acetylformoin were significantly enhanced in the presence of reductones such as ascorbic acid
  or methylene reductinic acid as well as the STRECKER-ACTIVE AMINO ACID proline."

  NEW STEP, family `Furanone_Amino_Acid_Reduction`, EXACT with NO token:
      hexose 1-deoxyosone + amino acid -> DMHF + Strecker aldehyde + CO2 + NH3 + H2O
      glycine:  C6H10O5 + C2H5NO2 -> C6H8O3 + CH2O  + CO2 + NH3 + H2O   C8H15NO7 both sides
      alanine:  C6H10O5 + C3H7NO2 -> C6H8O3 + C2H4O + CO2 + NH3 + H2O   C9H17NO7 both sides
  All six DMHF carbons still come from the hexose (the amino acid's carbon leaves as its own
  aldehyde and as CO2), so Wang & Ho's intact-C6 CAMOLA constraint is preserved — pinned in a
  test. The Strecker-aldehyde table was lifted to module level so this step and `_strecker_step`
  share one source of truth. The BARRIER IS UNCHANGED at 28.0, deliberately: it is the same
  value the step carried through `furanone_cyclisation`, so any movement in furaneol is
  attributable to the stoichiometry and not to a barrier change.

  THE MEASUREMENT THE BRIEF ASKED FOR (glucose + glycine, pH 5.5 / 150 C, max_generations=3,
  `_aminoketone_condensation` monkey-patched to return [], in-process, no repo file touched):
      Wave L1:  [baseline] 16 steps | DMHF steps 1 | H2 producers ['Aminoketone_Condensation']
                [pyrazine OFF] 14 steps | DMHF steps 0 | H2 producers []
      Wave P:   [baseline] 16 steps | DMHF steps 1 | H2 producers ['Aminoketone_Condensation']
                [pyrazine OFF] 15 steps | DMHF steps **1** | H2 producers []
  DMHF now SURVIVES. The `[HH]` token still exists and still has exactly one producer in an
  amine-only system; what changed is that furaneol no longer needs it. `_thiol_addition` (FFT)
  is still token-gated and is OUT OF SCOPE here. [P]

### (i) ITEM 7 — the verified Henry's-law constants

  NO CHANGE WAS NEEDED: the values are ALREADY ACTIVE and have been since 2026-08-26 (round 3,
  owner-approved; see the ledger entry of that date). `data/lit/henry_constants.yml` ships
  acrylamide `Kaw_25c: 6.0e-8` (`previous_value: 1.0e-9`) and HMF `5.0e-8`
  (`previous_value: 1.0e-10`), both above the `_NON_OBSERVABLE_KAW_THRESHOLD = 1.0e-8` gate in
  `src/recommend.py`, and both classify `observable`. No `literature_value_pending_threshold_
  decision` key survives anywhere.

  RE-MEASURED EMPIRICALLY TODAY, because "verified inert earlier" is a claim with a date on it.
  A throwaway harness (scratchpad, no repo file touched) scored the ENTIRE live benchmark panel
  twice in one process — once as shipped, once with both `Kaw_25c` values patched back to their
  `previous_value` placeholders in `recommend._HENRY_LOOKUP`:
      rows whose prediction differs: **NONE — 0 of 42 scored rows**
  So the activation remains numerically inert on the panel, and the earlier measurement holds.
  The reason is structural and worth writing down: acrylamide is produced by `src/safety.py`,
  which is not on the Henry-gated loop at all, and HMF's only appearances are in benchmarks
  where it is not a scored analyte. The partition layer as a whole is still inert (the `[D]`
  Henry workstream), so what these two values buy today is correctness of record, not
  correctness of prediction. Recorded so the next reader does not re-open it.

### (j) EVERY OTHER PANEL NUMBER THAT MOVED, old -> new

  BENCHMARK PANEL (`benchmark_summary.json`), max_ratio / MALE:
    cys_ribose_140C_Hofmann1998            2.2519 / 0.2192 -> **1.4110 / 0.0935**  (item 1)
    pea_isolate_uht_140C_Trikusuma2019     1.0000 / -      -> **2.2727 / 0.1189**  (item 4)
    resconi_2023_pbma_beef_identity        4.5607 / 0.6590 -> **4.6573 / 0.6681**  (species added
                                                                 to the shared volatile budget)
    thiamine_cys_glucose_120C_Bolton1994 772.8112 / 2.8881 -> **763.5881 / 2.8829**  (same cause)
    every other row bit-identical.
    NOTE ON THE TRIKUSUMA "BEFORE" MALE: the pre-Wave-P artifact carried 1.4562 dex next to
    max_ratio 1.0000, which is internally inconsistent (three rows all recovering exactly should
    give ~0). It is reported as read rather than reconstructed; the post-Wave-P 0.1189 is
    log10(2.2727)/3 and is self-consistent. Worth a look. [P]
    status_counts: scale-gap 7 -> **8**, pass 5 -> **4**, pass-no-ranking 2 UNCHANGED.
    evidence-role split 6/4/4 UNCHANGED. predictive passes **0/6 UNCHANGED**.
    fit-recovery passes 1/4 -> **0/4**. strict-ready **0/14 UNCHANGED**.
    TWO COUNTERS, AND THEY DISAGREE — found while re-pinning. `benchmark_summary.md` prints
    "0/6 + 2/4 + 4/4 = 6/14" because `src/presentation.py::_is_pass` also counts
    `pass-no-ranking`; the guard test and the README headline use strict `overall_status ==
    "pass"` and see 0/6 + 0/4 + 4/4 = 4/14. The divergence PREDATES this wave (7/14 vs 5/14
    after Wave O, where the ledger quoted only the strict number). Both are now pinned and the
    README says which is which. [P] reconcile the predicates.
  MC PANEL (`prediction_uncertainty.json`, n=200 seed 0): every COUNT unchanged — benchmark_count
    11, matched rows 35, coverage 29/35, honest_literature_coverage 1/3, not_evaluable 4,
    excluded_fitted_rows 2/2. Only interval widths moved: external_literature median CI
    0.8558 -> **0.8495** dex, fitted_row 2.3205 -> **2.2767**, internal_synthetic
    3.7657 -> **3.6929**.
  VALIDATION OVERVIEW: status_counts as above; inside_1_5x 3 UNCHANGED, outside_1_5x 6, worst
    quantitative point unchanged (CML 1203.68x).
  SNAPSHOTS: NOTES bumped v8 -> v9 with all six items and the honest movement. The four
    synthetic snapshots MOVED this time (unlike Wave O): the sulfur / pyrazine / furfural rows
    fall together by ~0.87x because the added species dilute the bounded volatile budget, and
    Nonanal falls by the oleate ratio. In the pea snapshot MFT rises to rank 1.
  PENTOSE >> HEXOSE: 3.39x -> **6.15x** (ribose 370.3 -> 686.83, glucose 109.3 -> 111.65).
    THIS IS NOT IMPROVED SUGAR DISCRIMINATION, and the test that pinned it said so in advance.
    Structural share re-measured in-process (thiol_addition_pentodiulose set equal to
    thiol_addition_hexose): **2.3118x**, up from 1.13x under Wave N. So 2.31x is mechanism and
    ~2.7x is the 3.30 kcal/mol gap between a FITTED barrier and an UNCONSTRAINED LEGACY FIT.
    More of the ordering rides on a fitted constant than before, even though that constant is
    better provenanced than the estimate it replaced. Both halves are in the README.
  PROJECTION REFIT (still NOT applied): fitted tau 5011.87 -> **10000** min, i.e. the budget
    scale the fit wants fell 2.51x -> **1.26x**; objective at fit 0.8620 -> 0.8817, at the
    shipped tau 0.8811 -> 0.8879.
    GENERATOR BUG FIXED WHILE HERE — the SECOND pass of the Wave M fix. Wave M derived the two
    residuals the decision prose quotes but left the ALLOCATION arithmetic hard-coded
    ("~1126 ppb", "542 ppb", "~51%", "~2800 ppb", "1.22x over", "~3.1x over", "2.25x MFT gap"),
    and every one had gone stale again. They are now computed from the live prediction on every
    run: tracked pool 1152.3 ppb over five species, furfural share 44%, MFT gap 1.41x under,
    FFT 1.09x over. The record's pointer to the stale `sulfur_barrier_refit_hofmann.md` was
    repointed at the Wave P record. This is Wave M's [P] item 3, closed one generator at a time.

### (k) TESTS RE-PINNED (each with a dated causal comment; none relaxed)

  1. `tests/unit/test_chemistry_soundness.py::test_pentose_reaches_norfuraneol_without_a_
     reduction_hexose_does_not` — the hexose probe was testing the `[HH]` GATE, not the
     chemistry, and the gate was the defect. Rewritten to probe the reductant the chemistry
     actually uses: pentose needs none, hexose+`[HH]`-token STILL emits nothing (asserted, so
     the gate cannot come back), hexose+amino-acid emits `Furanone_Amino_Acid_Reduction`.
     STRICTER than before: it now pins the absence of the gate as well as the presence of the
     route.
  2. `tests/unit/test_wave_h_2026_08.py::test_corrected_mft_route_barriers_are_estimates_not_
     fits` — `thiol_addition_pentodiulose` is no longer ESTIMATED/UNCONSTRAINED, so the old
     assertions would have been asserting a lie. Re-pinned to 26.35 with FITTED provenance, the
     verbatim conversion caveat, an explicit `!= 26.85`, and `deoxyosone_reduction` still
     ESTIMATED/UNCONSTRAINED.
  3. `tests/scientific/test_free_aa_quantitative_regression.py` — MFT band (1.85, 2.252, 2.80)
     -> (1.15, 1.411, 1.75); FFT (1.00, 1.219, 1.60) -> (1.00, 1.090, 1.43). Band widths held at
     the same RELATIVE span, so this is a re-centring and not a loosening. The header records
     that the contract now fails on ONE criterion (MALE 0.0935 > 0.09) rather than two, and that
     the improvement is FIT RECOVERY on a declared fit target.
  4. `tests/scientific/test_honest_headline_guards.py` x4 — fit-recovery passes 1 -> 0 with the
     full oleate/linoleate explanation; aggregate 5 -> 4 plus a NEW assertion that every
     remaining pass is internal-synthetic; a NEW pin on the lenient presentation-layer count
     (6/14) so the two predicates cannot drift apart unnoticed; MC median CI 0.8558 -> 0.8495
     with the companion widths; pentose ordering 3.39 -> 6.15 with the structural decomposition
     measured rather than asserted; and the hold-out guard gains a dated paragraph on the two
     nonanal points WITHOUT changing a single pinned number, because none of them moved.
  5. NEW `tests/scientific/test_wave_p_chemistry_2026_08.py` (20 tests) — every new family has
     an explicit barrier key and canonicalises to ITSELF (the fallthrough guard); none picks up
     a pH/aw correction from its name; the refit value and its caveat; the Hofmann rows; the
     C2+C3 lane's existence, its hexose reachability, and its ZERO contribution; norfuraneol's
     consumers and the isomer it must not emit; the oleate pool; the Trikusuma 2.2727x = 1/0.44
     identity; the fructose route; the pyrazine-independence measurement; the token-free DMHF
     step; and the intact-C6 CAMOLA constraint.

  NOT RE-PINNED, DELIBERATELY: `cys_ribose_140C_Hofmann1998`'s own contract (1.45x / 0.09 dex),
  the Pratap-Singh ratio tolerance, the Trikusuma observability factors, the 3.0x pentose/hexose
  floor, `max_fold_error`, the pre-widening 1/5, the 0/6 predictive headline, strict-ready 0/14.

### (l) [P] CARRIED FORWARD

  1. `src/recommend.py` selects the lowest-span path per product and does not sum parallel
     channels. Every compound claim rests on one route. Fixing it is a model-wide recalibration.
  2. `src/conditions.py` selects pH/aw corrections by SUBSTRING on the family name, so a
     family's NAME is a load-bearing calibration input. Wave P worked around it and pinned the
     workaround; it did not fix it.
  3. The C2+C3 lane does not express the low-moisture dependence Hofmann measured it under. Any
     gate would need an invented threshold or a naming trick.
  4. The oleate pool uses the linoleate rate constant; predicted nonanal stays biased HIGH by
     roughly an order of magnitude. Needs its own anchor.
  5. 2-mercapto-3-pentanone has no retrievable peer-reviewed odour threshold and is therefore
     not a scored sensory target.
  6. The two panel pass-counters (strict vs presentation-layer) disagree, 4/14 vs 6/14.
  7. The Trikusuma pre-Wave-P MALE (1.4562 next to max_ratio 1.0000) is internally inconsistent.
  8. `docs/validation/isotope_topology_evidence.md` mis-attributes 10.1021/jp0500900 to
     "Hearn & Smith (2007)"; CrossRef says Hung, Katrib & Martin (2005). Dossier not edited
     (do-not-touch list); the code uses the correct attribution.
  9. `_thiol_addition` (FFT) is still `[HH]`-pool-gated. Out of scope for item 6.
 10. The hexose retro-aldol channel that Yaylayan & Keyhani show gives glycolaldehyde directly
     from glucose is still absent; the C2+C3 lane reaches hexose systems by a longer route.

### (m) GATES + SUITE — and an honest note about the suite

  GATES, all three re-run LAST, after every code, artifact, test and doc edit:
    citation_gate  PASS — 81 files, 951 DOI-bearing fields, 316 unique DOIs (was 79/904/310
                   at Wave O). The six new DOIs passed with no formatting change; WAIVERS and
                   TEXT_SURFACE_WAIVERS are both still the empty tuple.
    holdout_guard  PASS — 3/3 invariants, re-run after the regeneration.
    fit_target_gate PASS — both checks. It now enumerates the new record
                   `sulfur_barrier_refit_pentodiulose.json: leverage=per_row_recovery,
                   1 target(s)` alongside the Wave H and Wave O records.

  FINAL SUITE: **1265 passed, 1 skipped, 2 xfailed, 0 FAILED** in 873.99 s
  (`tests/unit tests/scientific tests/integration tests/scripts`, documented conda path,
  exit code 0, zero FAILED/ERROR lines). Arithmetic check: 1265 + 1 + 2 = **1268**, which is
  the `--collect-only` count, and 1268 = Wave O's certified 1248 + the 20 new tests in
  `tests/scientific/test_wave_p_chemistry_2026_08.py`. The 1 skip and both xfails are the
  declared, strict-marked ones from Wave J2 (`xfail_strict = true` is on, so neither xfail can
  silently start passing).

  HOW THAT NUMBER WAS OBTAINED, because the route to it matters for trusting it. Three earlier
  attempts at this suite were abandoned, and the reason was the MACHINE, not the tree:
    * Attempt 1 was contaminated. The machine's data volume reached **100 % capacity**
      mid-run (133 MiB free at the worst point) and ENOSPC produced a cluster of spurious
      ERRORs at 17 % in subprocess-backed unit tests. Those do NOT reproduce once space was
      freed (clearing `__pycache__` and `.pytest_cache` recovered several GB), and the same
      region ran clean in every later attempt. Anyone reading a failing log from this period
      should check `df` before bisecting.
    * Attempt 2 reached **82 % with exactly ONE failure**, located by index against a
      `--collect-only` listing:
      `tests/scientific/test_matrix_only_benchmark.py::
       test_trikusuma_heated_pea_matrix_fit_is_recovered_to_within_1_05x`.
      That was a STALE PIN from item 4 — the nonanal recovery this wave deliberately breaks —
      and it is the one real failure the wave produced. Re-pinned two-sided with a dated
      causal comment; the file passes 7/7. All of `tests/unit` (the first 72 %) was green in
      that attempt.
    * Attempts 2 and 3 then stalled with free RAM at ~60 MB: pytest held ~99 % CPU while
      advancing a few tests per ten minutes. This is an environment condition and is not
      caused by anything in Wave P.
  The certified run above is the one that completed after the Trikusuma re-pin landed, on the
  final tree, with memory recovered. Its ZERO failures are what closes out the wave.

  TARGETED RUNS MADE ALONG THE WAY, all green, kept because they isolate the blast radius:
    - the 13-file Wave P blast radius (wave_p_chemistry, matrix_only_benchmark,
      matrix_assertions, matrix_headspace_ph_validation, honest_headline_guards,
      free_aa_quantitative_regression, pentose_hexose_sulfur_ordering, benchmark_summary,
      matrix_experiment_intake, lipid_oxidation_saturation, lipid_oxidation_guard,
      deep_research_benchmark_wiring, matrix_observability_refit_wave_o): **64 passed**
    - tests/unit/test_chemistry_soundness.py + tests/unit/test_wave_h_2026_08.py: **42 passed**
    - tests/unit/test_data_ingest.py + tests/scientific/test_matrix_experiment_intake.py:
      **5 passed**
    - tests/scientific/test_matrix_only_benchmark.py after the Trikusuma re-pin: **7 passed**

  NOT COMMITTED, NOT STASHED — handed to the orchestrator as instructed.

### (n) ENVIRONMENT FINDING (not a repo defect, but it invalidated a test run)

  `/System/Volumes/Data` hit **100 % capacity** during this wave (133 MiB free at the worst
  point). The observable symptom inside pytest was a cluster of ERRORs in subprocess-backed
  unit tests, which looks exactly like a code regression and is not one. Clearing
  `__pycache__` and `.pytest_cache` from the repo recovered several GB on its own. Anyone
  reading a failing CI log from this period should check `df` before bisecting.

## Wave P certification (2026-08-27, orchestrator)

The [P] "re-run the full suite on a machine with headroom before committing" is CLOSED:
full 4-directory suite under the final Wave P tree: 1265 passed, 1 skipped, 2 xfailed,
0 failed (873.99 s), run after freeing 6 GB of stale Docker images (the disk-full
condition that produced Wave P's spurious ERRORs and stalls). tests/integration and
tests/scripts are hereby certified; the exit-144 kills during the stall window were
Wave P's own process cleanup racing the orchestrator's runs, not test failures.

---

## Wave R — the Liu hold-out correction + the cleanup execution (2026-08-27)

Two parts. Part A is a data-honesty correction to reference values (the hold-out stayed out
of every fit; nothing was calibrated). Part B executes the Wave R0 cleanup inventory.

### PART A — THE LIU HOLD-OUT WAS BEING GRADED AGAINST NUMBERS FROM NOWHERE

Wave O's sharpest carried `[P]` — "Bi 2020 (1260 ppb) vs Liu 2023 (15-180 ppb band) disagree
by 24x on nominally the same system and no observability factor satisfies both" — is closed,
and the resolution is that ONE SIDE OF THE CONTRADICTION WAS NOT IN THE LITERATURE.

Source located, retrieved and READ IN FULL: Yaozheng Liu, "Flavor Chemistry of Pea Proteins",
MS thesis, North Carolina State University, **2021** (NC State Institutional Repository item
db647868-5ffe-4621-9f11-bbc4db357406); published as Liu, Cadwallader & Drake (2023), Food
Chemistry 406:134998, `10.1016/j.foodchem.2022.134998` (CrossRef-verified). Method: 24
commercial pea proteins rehydrated to 10% solids (w/w) in DI water; nine of them quantified
in Table 2.7 by HS-SPME-GC-MS/MS against five-point EXTERNAL standard curves with response
factors (R^2 > 0.95). Table 2.3 identifies the nine: samples 4, 12, 13, 14, 16 (ISOLATES) and
1, 2, 19, 22 (CONCENTRATES).

#### (a) PER-COMPOUND CORRECTIONS, OLD -> NEW, WITH THE SOURCE ROWS QUOTED

  **1. `liu_2023_ppi_hexanal_band` — hexanal. CORRECTED.**
    OLD: band 15.0-180.0 ppb, geometric mid **51.96**, OAV 3-40. Matched NOTHING; no hexanal
    number in the thesis is within two orders of magnitude of it.
    THESIS Table 2.7, row `Hexanal 2` (column order 1, 2, 4, 12, 13, 14, 16, 19, 22):
      `4318 +/- 481, 3360 +/- 274, 2445 +/- 212, 6052 +/- 776, 6383 +/- 409,
       11203 +/- 2216, 12181 +/- 4237, 52454 +/- 34133, 2533 +/- 804` ug/L
    NEW: band **2445-52454**, geometric mid **11320** (4 s.f.), span 21.45x, OAV 543-11656
    (Table 2.8, against the thesis's own 4.5 ug/L water threshold — which equals this repo's
    hexanal ODT, so the OAVs transfer without rescaling). The old band was 50-300x LOW.
    RECORDED, NOT SCORED: the isolate-only subset (4, 12, 13, 14, 16) is 2445-12181, mid 5457.
    The 52454 top of the scored band is sample 19, a CONCENTRATE and the lowest-protein
    sample of the panel (52.8 wt% protein, Table 2.3).

  **2. `liu_2023_ppi_nonanal_band` — nonanal. CORRECTED.**
    OLD: band 5.0-50.0 ppb, geometric mid **15.81**, OAV 5-50. Matched NOTHING.
    THESIS Table 2.7, row `Nonanal 2`:
      `0.797 +/- 0.55, 0.469 +/- 0.30, 0.188 +/- 0.040, 0.301 +/- 0.075, 0.652 +/- 0.11,
       0.861 +/- 0.11, 1.06 +/- 0.77, 1.15 +/- 0.66, 3.42 +/- 0.71` ug/L
    NEW: band **0.188-3.42**, geometric mid **0.8018**, span 18.19x. The old band was 6-266x
    HIGH. INDEPENDENT CORROBORATION IN THE SAME DOCUMENT: Table 2.8 reports nonanal OAV as
    `<1 <1 <1 <1 <1 1 1 1 3` against a 1 ug/L threshold — only consistent with sub-ppb-to-3-ppb.

  **3. `liu_2023_ppi_heptadienal_band` — (E,E)-2,4-heptadienal. RETIRED, NOT REPAIRED.**
    OLD: band 0.5-8.0 ppb, mid 2.0, OAV 7-114.
    THE COMPOUND IS ABSENT FROM THE SOURCE. Case-insensitive grep for `heptadien` over the
    full thesis text returns **ZERO hits** — not in Table 2.7 (quantitation), 2.4 (SPME-GC-O),
    2.5 (SAFE-GC-O) or 2.6 (AEDA/FD). The unsaturated aldehydes the thesis DOES quantify are
    (E,E)-2,4-nonadienal (2.05-19.1), (E,E)-2,4-decadienal (0.436-8.29) and (E)-2-decenal
    (3.65-143) ug/L. They are recorded in the correction note and DELIBERATELY NOT
    substituted: different compounds, different thresholds, and swapping one in would repeat
    the original error in a new form. Marked `value_status: no_verifiable_source`, removed
    from any scored target. It was never in the scored bundle, so no benchmark number moves.

  **4. `liu_2023_ppi_ibmp_band` — 3-isobutyl-2-methoxypyrazine. RETIRED, WRONG COMPOUND.**
    OLD: band 0.0-0.08 ppb, ODT 0.002, OAV floor 1.0.
    IBMP IS NOT IN THE SOURCE. The thesis states its pyrazine inventory explicitly: "The
    pyrazines identified in this study included 2-ethyl-6-methyl-pyrazine, IPMP,
    2-isoamyl-6-methylpyrazine, 2,5-dimethyl-3-(3-methylbutyl)-pyrazine,
    2,3-diethyl-5-methyl pyrazine and 2-sec-butyl-3-methoxy pyrazine." IBMP appears only in
    the literature review, quoting Murray et al. (1970/1976) on GREEN PEAS.
    NEAREST REAL ROW, RECORDED BUT NOT SUBSTITUTED: Table 2.7 row
    `2-Isopropyl-3-methoxypyrazine 2` = `57.0, 23.3, 18.2, 15.7, 21.9, 19.7, 29.0, 14.3,
    6.126` ug/L, i.e. an IPMP band of **6.126-57.0** (mid 18.69), Table 2.8 OAV 3063-28500
    against a 0.002 ug/L threshold. That is **76x-713x ABOVE the retired 0.08 ppb ceiling**:
    the record was not merely mis-labelled, it was three orders of magnitude wrong even for
    the compound it was nearest to. A properly-typed IPMP anchor is a separate owner decision.

  **5. `liu_2022_ppi_oav_anchors` (`data/lit/computational_priors.json`). CORRECTED.**
    Same underlying paper (the registry's own `shared_anchor_note` says so). OLD:
    `hexanal_oav_ppi: 28.0`, note "hexanal OAV=28 (dominant off-note), 1-octen-3-ol OAV=14",
    rank `[hexanal, 1-octen-3-ol, nonanal, pentanal]`. ALL WRONG. Table 2.8: hexanal OAV
    `960, 747, 543, 1345, 1418, 2490, 2707, 11656, 563` (19x-416x the retired 28.0);
    1-octen-3-ol `11, 7, 6, 4, 5, 3, 2, 4, 47` (never 14); 2-pentylfuran `89, 18, 13, 25, 40,
    38, 60, 4, 271` — which outranks 1-octen-3-ol in 8 of 9 samples and was MISSING FROM THE
    RANK LIST ENTIRELY. NEW rank `[hexanal, 2-pentylfuran, 1-octen-3-ol, nonanal, pentanal]`,
    `hexanal_oav_ppi` 960.0 with `hexanal_oav_ppi_range [543, 11656]`. NOT CONSUMED BY `src/`
    (only the archived one-shot ingest script that wrote it), so no prediction moves.

  **6. `liu_2023_ppi_offnote_baseline.key_values` (`data/lit/benchmark_intake_registry.json`).**
    A verbatim copy of the four defective bands. CORRECTED to the Table 2.7 ranges, plus
    `2_pentylfuran_ppb_range [25.5, 1628.0]` and `methional_ppb_range [15.5, 182.0]` (both
    verified from the same table), `ee_2_4_heptadienal_ppb_range -> null`,
    `methoxypyrazine_ppb_max 0.08 -> 57.0` (IPMP), `replicate_lot_count 6 -> 2`.
    FLAGGED NOT CORRECTED: `meaty_potential_multiplier: 0.12` is a repo construction the
    source does not report. `mft_detected: false` is consistent with the source (no
    2-methyl-3-furanthiol in any of its tables). **[P]**

  EVERY corrected record carries a `content_correction_note` with: the previous value, the
  statement that it matched nothing, the thesis table and row quoted verbatim, the UNITS
  BASIS (ug/L of the rehydrated 10%-solids slurry, read as ppb at rho ~ 1 kg/L — NOT a
  per-gram-of-powder or headspace basis), and the DI-WATER CALIBRATION CAVEAT (thesis 2.9:
  "pea protein solutions were replaced with 5 mL DI water and were spiked with known
  compounds at five different levels" — matrix binding is uncorrected, so **every value is a
  LOWER BOUND on total analyte**). `src/external_validation.py` gained a 4-line pass-through
  so `content_correction_note` travels from the payload anchor into the generated bundle;
  without it a regeneration would present the corrected number as if it had always been right.

#### (b) CITATION FIXED

  `"Liu, Y. (2023 thesis)"` -> `Liu, Y. (2021), "Flavor Chemistry of Pea Proteins", MS thesis,
  North Carolina State University; published as Liu, Cadwallader & Drake (2023), Food
  Chemistry 406:134998`. The old string was doubly wrong: the thesis is **2021** (the 2023
  belongs to the derived paper). `doi` added to all four payload anchors and `source_doi` to
  the bundle: `10.1016/j.foodchem.2022.134998`. This REVERSES the Wave I note that said the
  DOI was "deliberately NOT used here" — that reasoning confused CITING the paper with
  PROMOTING the hold-out. Evidence class `external_validation_only` is unchanged and the
  bundle remains outside every fit. The record ids keep their `liu_2023_` prefix, which is
  correct for the published version and is referenced by four tests and two backlog files.

  A STATUS-SEMANTICS DISTINCTION, made deliberately: the two retired anchors carry
  `source_status: verified_primary_source_read` + `value_status: no_verifiable_source`, NOT
  `source_status: no_verifiable_source`. The CITATION is verified — the source was identified,
  retrieved and read, and its DOI is real. What has no verifiable source is the VALUE. Writing
  `source_status: no_verifiable_source` beside a real DOI is internally incoherent and the
  citation gate's status-coherence rule catches exactly that. The gate was NOT relaxed; the
  data was made accurate. `value_provenance`, `pipeline_role` and `benchmark_role` all still
  read `no_verifiable_source` / `retired_unverifiable`, so a grep still finds these records.

#### (c) THE RANKING CONTRACT, RE-DERIVED

  Auto-derived from the corrected panel (no explicit `comparison_contract` is supplied):
  adverse markers `[hexanal, nonanal]`, hexanal expected_rank 1, nonanal 2, both direction
  `lower`. UNCHANGED IN SHAPE, and now for a real reason: on the corrected panel hexanal
  (11320) leads nonanal (0.8018) by 14000x, where the retired bands had it leading by only
  3.3x. The contract's `citation_provenance` now carries the corrected citation.

#### (d) THE HOLD-OUT, OLD -> NEW. IT GOT WORSE. NO PREDICTION MOVED.

  `results/validation/external_validation_report.json`, n=200 seed 0, hold-out never fitted.

  | point                         | meas old | meas new | p50      | fold old  | fold new  |
  | ----------------------------- | -------- | -------- | -------- | --------- | --------- |
  | bi_2020_raw_pea / hexanal     |   1260.0 |   1260.0 |   1013.1 |    1.2437 |    1.2437 |
  | bi_2020_roasted_pea / hexanal |    324.0 |    324.0 | 801700.4 | 2474.3839 | 2474.3839 |
  | li_2026_hme / 1-hexanol       |    20.04 |    20.04 |  22394.3 | 1117.4799 | 1117.4799 |
  | li_2026_hme / 2-pentylfuran   |   5625.8 |   5625.8 |  11038.6 |    1.9621 |    1.9621 |
  | li_2026_hme / hexanal         |    605.6 |    605.6 |  56409.9 |   93.1471 |   93.1471 |
  | li_2026_hme / nonanal         |    72.66 |    72.66 |   8596.4 |  118.3093 |  118.3093 |
  | **liu_2023_ppi / hexanal**    |    51.96 | **11320**|   1013.1 |   19.4972 | **11.1739** |
  | **liu_2023_ppi / nonanal**    |    15.81 |**0.8018**|    75.55 |    4.7785 | **94.2234** |

  * median |log10| 1.6296 -> **1.9717** dex; median fold 42.6159x -> **93.6837x**
  * shipped-sigma coverage 4/8 -> **3/8** (0.500 -> 0.375). The hexanal point stays INSIDE
    the interval; the nonanal point leaves it.
  * pre-widening (ln-sigma 2.0) coverage **3/8 UNCHANGED** — the two priors now AGREE,
    because the last point separating them left both intervals.
  * `genuine_extrapolation` **1/5 UNCHANGED**; `max_fold_error` **2474.3839x UNCHANGED**;
    `matched_compound_count` 8 UNCHANGED. `in_panel_rescoring` hits 3/3 -> 2/3.
  * `external_failing_compounds`: nonanal mean |log10| 1.3762 -> **2.0236**; hexanal 1.6868
    -> 1.6264 (IMPROVED); 1-hexanol 3.0482 unchanged.
  * SIX POINTS ARE BYTE-IDENTICAL. Only the two Liu MEASURED values moved. **Nothing was
    fitted, no constant was touched, no prior was changed.**

  THE NONANAL ROW, REPORTED PLAINLY. The model predicts **75.55 ppb against a real inter-lot
  band whose TOP is 3.42 ppb** — 94.22x over, and it is now the sharpest lipid-lane
  over-prediction the repo has against a DIRECTLY-QUANTIFIED reference (the larger folds on
  record — 2474x, 1117x, 118x — are all against reported point values in other lanes, and the
  roasted-pea 2474x is an OAV-derived number). WAVE P'S OLEATE FIX IS THE PARTIAL MITIGATION
  ALREADY LANDED: against the corrected reference the same fix reads **214.1x -> 94.22x**
  (p50 171.70 -> 75.55), i.e. it removed more than half the error in log space and was still
  not close. And the DI-water caveat runs the wrong way for us: Liu's 0.188-3.42 is a LOWER
  bound on total nonanal, so the true over-prediction is if anything larger than 94x.

  WAVE O IS VINDICATED BY THE CORRECTED TARGET. Wave O accepted a headline regression on the
  grounds of a 24x literature contradiction. There was no contradiction: 51.96 was a midpoint
  of a band that appears in no source. Against Table 2.7 the verified 1138 ppb Pratap-Singh
  anchor sits just UNDER Liu's lowest lot (2445) rather than 6.3x above her band, Bi's 1260
  and Liu's 2445-52454 are consistent for different matrices, and this same Liu hexanal point
  IMPROVED 19.50x -> 11.17x once graded against the real number. The refit was right and was
  being marked by a broken key.

#### (e) A REPRODUCIBILITY DEFECT FOUND WHILE REGENERATING (fixed)

  The committed `external_validation_li_2026_*.{yaml,json}` did not reproduce from their own
  generator. Two causes, both fixed:
  1. `_measurement_row`'s Wave I 4-significant-figure rounding was applied to EVERY number,
     including `reported_point_value` rows. It was written to strip precision WE invented
     (sqrt(15*180) = 51.96152422706632); applied to a measurement it strips precision the
     SOURCE has. Li 2026 2-pentylfuran was being silently rewritten 5625.8 -> 5626.0 on every
     regeneration, and Wave R's corrected Liu band max would have been written 52454 -> 52450.
     Rounding now applies ONLY to constructed values (band midpoints, OAV-derived points);
     never to reported values, and never to band endpoints, which are per-lot measurements
     quoted from the source table. Dated causal comment in the docstring.
  2. Wave K/M's Li 2026 correction notes lived inline in `value_provenance_note` in the
     GENERATED files but as a structured `content_correction_note` in the payload, which the
     generator dropped. The new pass-through surfaces them properly. No number changed.

#### (f) TESTS RE-PINNED (each with a dated causal comment; none relaxed)

  1. `tests/unit/test_audit_remediation_carried_2026_08.py::test_holdout_bundle_carries_typed_
     identifier_and_keeps_its_evidence_class` — `assert "source_doi" not in bench` INVERTED to
     `assert bench["source_doi"] == "10.1016/j.foodchem.2022.134998"`, i.e. STRICTER: the DOI
     is now required to be that specific one, so the field cannot drift back to empty or to
     another paper. Two new assertions require the identifier to name the thesis and the year
     2021. The evidence-class invariant is untouched — that is the point of the test.
  2. `tests/scientific/test_honest_headline_guards.py::test_holdout_scores_1_of_5_...` —
     `median_accuracy_fold` 42.62 -> **93.68**, `ci_coverage_hits` 4 -> **3**,
     `ci_coverage_rate` 0.5 -> **0.375**, each with the full two-point movement table and the
     statement that no prediction moved. `max_fold_error` (2474.4), the pre-widening 1/5, the
     5/3 kind split and `matched_compound_count` 8 are DELIBERATELY NOT re-pinned — they did
     not move, and their stability across a reference correction is the guard's whole value.
     The docstring gains a block recording that its own quoted "Liu 4.78x (IMPROVED)" was
     computed against a fabricated reference.

  NOT RE-PINNED, DELIBERATELY: `test_slr_reference_payloads.py:131` and
  `test_literature_learning_loop.py:77` still map the backlog entry to
  `liu_2023_ppi_ibmp_band`. The RECORD still exists (retired, not deleted), so the mapping is
  still true; deleting the id would delete the retirement record. `deep_research_backlog.json`
  keeps its raw `"Liu, Y. (2023 thesis)"` citation string with its existing repair record —
  that is the historical raw text, and `test_runtime_first_registry_landing.py:119` reads it.

#### (g) DOC SYNC (numbers match the regenerated artifacts exactly)

  README.md: hold-out surface row 42.62x -> 93.68x; the "median hold-out fold error is 42.62x"
  sentence -> 93.68x; a new leading Wave R block (old-vs-new bands, per-point folds, the
  retired compounds, the DI-water caveat, the explicit "no prediction moved"); the Wave O
  block's 24x-contradiction paragraph marked SUPERSEDED with the corrected reading; the Wave P
  block's "all unchanged" annotated with what those numbers are now and with the 214x -> 94x
  reading of the same fix. The two historical "15.31x -> 42.62x" statements are left as
  history, because they are.
  AUDIT.md: headline hold-out row; three rows in "What Round 3 cost"; a new "Wave R" section
  under Round 3 with the five-row old-vs-new table, the price table, and the vindication note;
  the Wave O 24x paragraph marked superseded in place.
  docs/reference/VALIDATION_CONTRACT.md section 3E: "Current numbers" re-headed to the Wave R
  regeneration, 42.62x -> 93.68x, 4/8 -> 3/8 (and the note that both priors now agree), a new
  Wave R paragraph, and supersession notes on the Wave O and Wave P paragraphs.

### PART B — CLEANUP EXECUTED (against `wave_r0_cleanup_inventory.md`)

  **TIER 1 — DONE.**
  * `src/Maillard.egg-info/` — 5 tracked files untracked and deleted (332 lines). Its
    `top_level.txt` still advertised the five modules the audit deleted, i.e. the repo was
    shipping a manifest asserting the dead-code island existed. `*.egg-info/` added to
    `.gitignore` with a dated causal comment.
  * `tasks/.todo.md.bak` — deleted (1162 lines). `*.bak` added to `.gitignore`.
  * `src/validation/__init__.py` — deleted; the empty orphan package is gone. Every
    `src.validation` hit in the tree is the MODULE `src/validation_contract.py`.
  * **All 14 dead top-level functions deleted (~186 lines), after RE-VERIFYING all 14 —
    including the two the inventory deferred because Wave P was editing those files.** Each
    name returns exactly one hit (its own `def`) across `*.py/md/json/yaml/txt/cfg/toml`,
    excluding `scratch/`, `.git/`, the ledger and AUDIT.md:
    `load_report_result` (reporting), `build_multi_run_family_lane_sensitivity_payload`
    (family_lane_sensitivity), `render_family_payload_coverage_markdown`
    (literature_family_registry), `export_matrix_target_panel` (matrix_targets),
    `_publication_year_proxy` (external_validation), `render_domain_warnings_markdown` +
    `render_domain_warnings_cli` (presentation), `iter_family_prior_entries`
    (matrix_prior_registry), `get_process_state_calibration_payload` (matrix_correction),
    `build_refinement_governance_artifact` (selective_refinement_governance),
    `_parse_xtb_energy` (xtb_path_quality), `load_optional_json_mapping` (artifact_io),
    **`_predicted_order_lookup` (benchmark_validation)** and **`_species_from_pool`
    (smirks_engine)** — the last two re-checked post-Wave-P and still zero-reference.
    ONE import became orphaned and was removed: `Iterable` in `family_lane_sensitivity.py`
    (verified against HEAD: it was used only by the deleted function). Every other unused
    import in those 13 files pre-dates this wave and was left alone.
    `_publication_year_proxy`'s docstring claimed it was "kept only so importers do not
    break"; there were none. Its retirement rationale is preserved as a comment at the call
    site that used to invoke it — a retired date-fabricator sitting in the module is a loaded
    gun, and the comment is the part worth keeping.

  **TIER 2 — DONE.** `scripts/get_details_fgh.py` + `get_details_fgh_output.txt` moved to
  `scripts/ingest/archive/` with `git mv` (history follows). `scripts/ingest/archive/README.md`
  gains a "Deep-research detail extraction" catalogue entry and its three false claims are
  corrected in place, each marked with the date and the reason: (i) "Most of these files are
  not tracked in git … nothing here has been `git add`ed" — false, `git ls-files` returns all
  21 scripts plus the README; the corrected text keeps the true half (only the two
  `scratch_*.py` files have their pre-move HISTORY); (ii) "`scripts/sync_backlog.py` —
  untracked" — false, it is tracked; (iii) the `get_details_fgh` bullet's "untracked one-shot
  … a candidate for this archive" — both halves settled. The catalogue entry also records a
  path caveat the inventory missed and that is WORSE than the `.parents[N]` problem the README
  already documents: this script resolves `ROOT = Path(".")`, i.e. relative to the CWD, not to
  its own location, so it only ever worked from the repo root and archiving changes nothing.

  **TIER 3.1 — SKIPPED, AND THE INVENTORY'S CASE FOR IT DOES NOT HOLD.** The instruction was
  to consolidate only if the churn is genuinely zero. Call-site churn would indeed be zero
  (private names retained, bodies delegate), but the recommendation is wrong on two counts
  that only show up on reading the actual bodies:
  1. **The shared helpers ALREADY EXIST and are simply unused.** `src/artifact_io.py` already
     defines `repo_root()` (identical body) and `load_json_mapping()`. There is nothing to
     add. Twelve modules ignore them.
  2. **`_load_json` is not one helper, it is TWO with different failure semantics.** Three
     modules (`chemistry_benchmark_validator`, `computational_gap_refinement`,
     `literature_family_registry`) use a 3-line form with **no existence check** — a missing
     file RAISES. Eight modules use a 4/5-line form that **returns `{}`** on a missing file.
     Merging them would silently convert a raise into an empty dict in three modules. In this
     repository an empty dict is a silent-degradation path, not an error path; that is a
     behaviour change wearing a refactor's clothes, and it is precisely the class of defect
     this audit exists to find. The `{}`-returning variant was ALREADY in `artifact_io` as
     `load_optional_json_mapping` — dead, and deleted in TIER 1 of this same wave. The
     shared-helper approach was already tried here and abandoned.
  3. The claimed saving is also arithmetically wrong for `_repo_root`: the body is 2 lines and
     a delegating body is 2 lines, so 12 delegations save **0 lines and add 12 imports** (net
     +12). The inventory's "-24 lines" for `_repo_root` and "-45 lines" overall do not exist.
  **[P] for the owner:** the real finding is that `artifact_io.repo_root` / `load_json_mapping`
  have 12 and 11 uncoordinated shadows. Fixing that is a deliberate migration with a semantics
  decision (which `_load_json` behaviour is correct where), not a mechanical dedup.

  **TIER 4 — DONE (D1-D11).**
  * D1 `.github/workflows/ci.yml` — the 5-line header describing `tests/benchmarks/` and
    `_lane_policy.py` as live-but-unrun rewritten to say the lane was DELETED (Wave J2, one
    day before the comment outlived it), what survives (the loader + its tracked `data/qm/`
    fixtures), and that its only consumer is its own parse test.
  * D2 `VALIDATION_CONTRACT.md` section 7 — the present-tense "`tests/benchmarks/` is
    intentionally skip-heavy today" rewritten: the directory does not exist, nothing skips
    there because nothing runs there.
  * D3 same file — `scripts/generate_benchmark_summary.py` -> `scripts/generators/…`.
  * D4 same file — the copy-pasteable `targets … cys_glucose_150C_Farmer1999.json competing`
    line REMOVED, not repointed. The file was deleted as fabricated and there is no equivalent
    competing-target benchmark; inventing one to keep an example alive is the failure mode the
    document exists to prevent. The removal is annotated in place.
  * D5 `docs/architecture.md` — `src/skala_refiner.py` removed from the module list; it does
    not exist and never appears in `git ls-files` (documented before it was written, then
    never written).
  * D6 `docs/protocols/pea_matrix_meaty_benchmark.md` — the "as specified in
    `docs/EXTERNAL_MATRIX_BENCHMARK_UNLOCK_REPORT.md`" pointer removed; that document is
    absent from the repo and from git history.
  * D7 `scripts/prepare_xtb_dirs.py` + `scripts/run_tier2_dft.py` — both user-facing "Did you
    run …?" error strings repointed to `scripts/generators/generate_mapped_geometries.py`.
    These are what a stuck user is told to run, which is why they matter more than a docstring.
  * D8 `scripts/generators/ingest_dft_c4_c5_results.py` — "Run scripts/run_dft_c4_c5.py first"
    replaced by a statement of the actual gap: no such runner exists in the working tree or
    anywhere in git history, so the inputs must be produced by hand.
  * D9 `generate_mace_training_data.py` + `generate_rmg_inputs.py` — self-referential
    pre-move paths in their own docstrings corrected to `scripts/generators/…`.
  * D10 `data/qm/phase33_barrier_benchmarks.json` — **APPENDED, NOT OVERWRITTEN**, as a new
    `runtime_coupling_correction_2026_08_27` key beside the original. The original is an audit
    record of what was true when it was checked. The correction makes the finding STRONGER:
    with `_lane_policy.py` gone, NOTHING gates on this file's existence, so the eighteen
    unverified numbers are fully inert.
  * D11 `VALIDATION_CONTRACT.md` — the `thermodynamic_gating_audit.{md,json}` sentence
    annotated: both artefacts are gitignored and absent; the command regenerates them.
  * D12 (README `results/first_run/report.md`) — LEFT ALONE. It is a generated-output example
    path, correct as an example, and the inventory graded it low severity.

  **NOT TOUCHED, per the inventory's `keep` grades and the brief:** `direct_sulfur_bonus`
  (3.3), `depth_bias_strength` (3.4), `scripts/calibrate_barriers.py`'s dangling-path tripwire
  (3.5), `sulfur_barrier_refit_hofmann.{json,md}` (3.6), `src/authority_benchmark_data.py`
  (3.7), the `_mft_pathway` hexose lump and the `Thiol_Addition_Norfuraneol` retirement
  guards. `scratch/`, `.claude/`, `data/geometries/` and
  `docs/validation/isotope_topology_evidence.md` untouched.

### GATES

  citation_gate **PASS** — 81 files, 959 DOI-bearing fields, 316 unique DOIs (was 79 / 904 /
  310 at Wave O; the growth is Wave P's additions plus Wave R's four new `doi` fields on the
  Liu anchors — the DOI itself was already in the repo, so unique DOIs moved with Wave P, not
  with this wave). WAIVERS still empty.
  holdout_guard **PASS** — 3/3 invariants, re-run after the correction and after regeneration.
  fit_target_gate **PASS** — both checks.

### ARTIFACTS REGENERATED

  `generate_external_validation_payloads.py` (4 protocol YAMLs + 4 benchmark JSONs),
  `generate_external_validation_report.py` (`external_validation_report.{json,md}` +
  `external_failing_compounds.json`), `generate_external_validation_inventory.py`
  (59 candidates / 8 executable / 50 narrative / 1 redundant, unchanged),
  `generate_literature_backlog.py` and `generate_literature_learning_loop.py` (both read
  `flavor_reference_payloads.json`; outputs gitignored, no tracked change). The calibration
  panel, `benchmark_summary`, `prediction_uncertainty` and the matrix sigma derivation are
  STRUCTURALLY unaffected: the hold-out is excluded from `get_benchmark_files()` and no
  constant moved.

### FINAL SUITE (Wave R)

  **1265 passed, 1 skipped, 2 xfailed, 0 FAILED** in 1028.55 s
  (`tests/unit tests/scientific tests/integration tests/scripts`, documented conda path,
  `-p no:randomly`, exit code 0, zero FAILED/ERROR lines). 1265 + 1 + 2 = **1268**, the same
  collected count Wave P certified — Wave R added no tests and removed none; it re-pinned two
  existing assertions (one of them INVERTED to be stricter). The 1 skip and both xfails are
  the declared, strict-marked ones from Wave J2.

  HOW IT WAS OBTAINED, because one run was discarded. A first pass was started before Part B
  was finished and was **killed by this wave at 62%** — deliberately, and only because it was
  this wave's own process. Eleven `src/` files were rewritten under it mid-run (the
  whitespace-minimisation pass described below), which makes any result from it meaningless:
  it stalled at 62% precisely where the edits landed. It was killed, the tree was frozen, and
  the suite above was run once, start to finish, against the final tree with nothing else
  touching the repository. No pytest process this wave did not start was signalled.

  A WHITESPACE NOTE on the dead-function deletions. The first deletion pass collapsed runs of
  3-6 blank lines to 2 wherever it found them, which produced unrelated hunks in
  `benchmark_validation.py` (18 whitespace lines), `smirks_engine.py` (7) and others, and
  stripped a trailing-newline-less EOF in `family_lane_sensitivity.py`. That churn was
  REVERTED and the deletions redone so each of the eleven pure-deletion files carries exactly
  ONE hunk: the function and the blank lines that separated it from its neighbour, with the
  original spacing preserved. A cleanup wave should not smuggle in a reformat.

  NOT COMMITTED, NOT STASHED — handed to the orchestrator as instructed. 48 tracked paths
  changed (+1139 / -2046 lines); the three untracked `data/geometries/` directories that were
  present at entry are untouched, as are `scratch/`, `.claude/` and
  `docs/validation/isotope_topology_evidence.md`.

### [P] CARRIED FORWARD FROM WAVE R

  1. **The nonanal lane is 94x over against a directly-quantified reference** (75.5 ppb vs a
     0.188-3.42 band whose top is 3.42). Wave P's oleate-substrate fix removed more than half
     the log error and stopped there. The remaining bias is the model treating oleate as
     being as oxidisable as linoleate — see `src.lipid_oxidation.MARKER_HYDROPEROXIDE_POOL`.
     This is now the highest-value lipid-lane target, ahead of the 1-hexanol factors.
  2. **A properly-typed IPMP anchor.** The source quantifies 2-isopropyl-3-methoxypyrazine at
     6.126-57.0 ug/L with OAV 3063-28500 — a real, verified, high-OAV pea off-note the repo
     currently has NO anchor for, because the slot it would occupy was filled by a fabricated
     IBMP band. Adding it is an owner decision (new record, not a repair).
  3. **The Liu bundle's proxied pH 6.0 is not the source pH** (Table 2.3: 6.3-7.3, mean ~6.8
     across the nine quantified slurries). Left unchanged because editing an executable
     condition changes the prediction and Wave R corrected reference data only.
  4. **`meaty_potential_multiplier: 0.12`** on `liu_2023_ppi_offnote_baseline` is a repo
     construction the source does not report, and unlike the bands it IS consumed
     (`src/matrix_correction.py` reads the field name from the source-profile registry).
     Not corrected because there is nothing in the source to correct it to.
  5. **`artifact_io.repo_root` / `load_json_mapping` have 12 and 11 uncoordinated shadows in
     `src/`**, and the `_load_json` shadows split into two groups with DIFFERENT failure
     semantics (raise vs return `{}`). Consolidation is a deliberate migration with a
     semantics decision per call site, not the mechanical dedup the cleanup inventory
     proposed. See TIER 3.1 above.
  6. **The Phase 3 authority lane is now consumerless.** With `tests/benchmarks/_lane_policy.py`
     deleted, `src/authority_benchmark_data.py`'s only importer is its own parse test, and its
     eighteen unverified `data/qm/phase33_*` numbers are inert. Revive or retire — an owner
     decision that Wave R only documented (D1, D2, D10).
  7. **Wave I's 4-significant-figure rounding was silently altering measured values** on every
     regeneration, and the committed `li_2026` artifacts had drifted away from their own
     generator as a result. Fixed here for the hold-out generator. Whether the same
     round-everything pattern exists in other generators was NOT audited. [P]

---

## Wave S2 — directional (ordinal) accuracy panel + three hard retrievals (2026-08-27)

Read-only on `src/` and `data/`. New files: `docs/validation/directional_claims_panel.yml`
(52 claims, 14 sources) and `docs/validation/directional_accuracy_report.md`. The scoring
runner was deliberately left in the session scratchpad, NOT committed — an in-repo
directional scorer becomes an optimisation target the moment it enters CI, and this panel's
value is that the model has never seen it. Nothing was tuned.

### Why: the repo has never measured the thing it is actually used for

Every validation artifact here measures ABSOLUTE accuracy, and reports it honestly as bad.
But the model's real use is ORDINAL — which sugar, which isolate, which way to move a knob.
Nobody had scored that. This wave built the evidence base and ran it.

### Headline: 18/29 (62%) on strictly independent claims

| Bucket | Agree / evaluable | Unevaluable |
|---|---|---|
| **Strictly independent (headline)** | **18 / 29 (62%)** | 8 |
| + system-overlap (ratios the fit cancels out of) | 6 / 6 (100%) | 0 |
| Screening total | 24 / 35 (69%) | 8 |
| Fit-adjacent (per-row fit targets, **excluded**) | 9 / 9 (100%) | 0 |

That gradient — 100% on fitted rows, 100% on system-overlap, 62% out of sample — is the
finding. Most rows are binary, so chance is ~50%: the directional edge is real but modest.

### The failures are NOT spread evenly. Two knobs carry almost all of them.

| Category | Score | | Category | Score |
|---|---|---|---|---|
| sugar_identity | **8/8** | | **ph** | **2/7** |
| temperature | **6/8** | | **moisture_aw** | **0/3** |
| additive_cysteine | 3/4 | | time | 2/2 |
| lipid_lane | 2/2 | | matrix_identity | 1/1 |

**Excluding pH and aw the headline is 16/19 (84%). Restricted to pH and aw it is 2/10 (20%)
— worse than a coin.** The model is a usable screening tool over sugar identity, precursor
loading, temperature and time, and an actively misleading one over pH and moisture.

### THREE CODE FINDINGS, verifiable by inspection (not fixed — `src/` not owned this wave)

1. **`ReactionConditions.get_ph_multiplier` is dead on the prediction path.** It encodes the
   documented physics (pyrazines alkaline-favoured, 0.2→8.2 sigmoid; furans/thiols
   acid-favoured, 4.9x boost) and is unit-tested for the correct signs in
   `tests/unit/test_conditions.py`. `grep -rn get_ph_multiplier src/` finds it ONLY in
   `kinetics.py`, `pathway_ranker.py`, `cantera_export.py`.
   `benchmark_validation._run_benchmark_recommendation` calls `conditions.get_rate_constant()`,
   which applies `_ionization_correction` + `_water_activity_correction` and **not**
   `get_ph_multiplier`. Written, tested, and never executed where it matters. **[P]**
2. **The pyrazine branch of `_ionization_correction` (`conditions.py:310`) is unreachable.**
   It keys on the substring `"pyrazine"` in the reaction-family name. The families the
   SMIRKS engine actually emits are `Aminoketone_Condensation`, `Enolisation_2_3`,
   `Enolisation_2_3_Amadori`, `Retro_Aldol_Fragmentation`, … — none contains "pyrazine", so
   the branch returns 1.0 at every pH. Net effect: the model moves 2,5-dimethylpyrazine the
   WRONG WAY with pH in every system tested, against two independent direct measurements
   (Laemont & Barringer 2023 measure 26.6 → 37.4 → 68.2 ppb over pH 4→7→9; model gives
   99.5 → 16.8 → 14.9). **[P]**
3. **`_water_activity_correction` is correctly shaped and effectively inert.** It implements
   Labuza with a real peak at aw 0.65 and a 5x range — but is applied to only 2 of the 7
   emitted families (`Amadori_Rearrangement`, `Strecker_Degradation`), neither on the
   furfural/HMF track. At the observable the residual is a 4.9% wiggle (HMF 570 → 598 → 594
   over aw 0.25 → 0.65 → 0.95), below any decision-relevant threshold. This is a ROUTING
   fix, not a modelling one. **[P]**

Also: **the model can only make acrylamide go up with temperature** (TEMP-01, ACR-02). Ma
2024 measures 130 > 150 > 170 °C; a review places the acrylamide maximum at 160–180 °C. The
model has no acrylamide destruction term. For a tool whose pitch includes process-safety
screening, this is the most consequential single miss in the panel.

### Eight measured scope limits (each probed against the running model, not assumed)

* **Ribose and xylose are indistinguishable** — identical predictions to 4+ s.f. on every
  observable. Ni 2021 measures xylose giving 2.3x the FFT of ribose while ribose gives 1.8x
  the MFT of xylose. Any pentose recommendation is a coin flip on this axis.
* **No α-dicarbonyls, no CML, no CEL** anywhere in the output — despite
  `cml_cel_commercial_pbma_Foods2023.json` existing as a benchmark.
* **No methionine Strecker branch** — methionine + glucose at 140 °C returns ZERO named
  volatiles. No methional, though it is one of the largest responses in Trikusuma's table.
* **`MATRIX_BENCHMARK_PROFILES` has two entries.** Brown rice, which carries 14–20x the
  hexanal of either legume isolate in Pratap-Singh Table 1, cannot be ranked at all.
* **The matrix lane and free lane share no observable** — a protein isolate emits only the
  four lipid markers, so "what happens to pyrazines in a pea beverage" is unanswerable.
* **No off-note can go DOWN with heat.** Trikusuma measures beany 3.3 → 1.5 while
  oxidized/painty rises 1.0 → 2.7 under the same UHT step; all four matrix markers rise
  monotonically with thermal load.
* **2-Pentylfuran has no Maillard route**, only a lipid one, so the measured hexose > pentose
  ordering for it cannot be represented.

### MISSION 2 — three hard retrievals

**2a. Trikusuma — CONFIRMED.** `pea_isolate_uht_140C_Trikusuma2019.json`'s 782 / 163 / 24 ppb
are real. Source: **Trikusuma, M., MS thesis, Ohio State University 2018**, OhioLINK ETD
accession `osu1531495328317918`, open access — the direct precursor of the paywalled paper.
Table 6, p.35: hexanal 781.72 ± 58.59, 2-pentylfuran 163.16 ± 15.06, nonanal 23.98 ± 0.80
µg/L, "Processed" column. The ±8/9/4% in the benchmark match the SD-derived RSDs (7.5/9.2/3.3%).
Conditions match (140 °C, 6 s hold, pH 7.16 processed). **Two curation notes:** the canonical
citation year is **2020** (Food Chem 312:126082), not 2019; and the values are the Processed
column, not Control (331/59/8) or Aged (683/197/22) — the `heated_matrix` state label is
load-bearing and deserves an explicit note in the file. [O]

**2b. Hofmann & Schieberle 1998 — STILL OPEN.** Every route exhausted; body closed
everywhere (ACS 403, Unpaywall `oa_locations: []`, OpenAlex/S2 CLOSED, 5 OA citing papers
pulled in full text — none reproduces the yields, Google Books zero volumes, Cerny 2015
cites it 5x and quotes no number). Only quantitative text obtainable is the abstract, whose
mol % figures (1.4 / 0.05) belong to a DRY 180 °C intermediate system, not the benchmark's.
**Arithmetic, written out:** MFT/FFT MW 114.17; on the assumed 10 mM basis the repo's ppb
imply **MFT 0.0300 mol %** and **FFT 0.0175 mol %** (342e-6/114.17/0.010×100). Forward, the
abstract's 1.4 mol % would give 15 984 ppb (46.7x) and 0.05 mol % would give 571 ppb (2.85x).
Not confirmed, not refuted. **Three things for the owner:**
  1. The implied yields are chemically plausible but trace to nothing.
  2. **The conditions appear misattributed.** Cerny 2015 verbatim: Hofmann & Schieberle's
     aqueous ribose/cysteine model was **145 °C / 20 min at 1:3 cysteine:ribose**; the
     benchmark's **140 °C / 30 min / equimolar** is the **Mottram & Nobrega 2002**
     (`10.1021/jf0200826`) protocol — and they reported headspace GC-MS, not solution ppb,
     so they are not the source of 342/200 either.
  3. **The 10 mM basis is itself unverified**, so the conversion carries a free multiplicative
     parameter — sitting underneath the panel's TIGHTEST contract (1.45x / 0.09 dex).
  Two clean resolutions only: ILL/institutional access to the body, or re-express the target
  in the paper's native mol % and drop the ppb. `data/` not edited (concurrent agent owns it). [O]

**2c. Brands & van Boekel 2001 — RECOVERED. The "unreadable scan" finding was WRONG.**
`data/lit/timeseries/brands_sugar_casein_120C_pH68.yml` asserts the multiresponse data
"exist ONLY as figures in a 1-bit 200 dpi scan… reading points off those figures would be
eyeball estimation and was NOT attempted." The scan IS bitonal CCITT at 200 dpi, but
rendered at 400 dpi the markers are 20–25 px across and unambiguous — verified directly by
rendering the page and inspecting the annotated overlay. **Also: DOI `10.1021/jf001430b` is
"Reactions of Monosaccharides during Heating of Sugar−Casein Systems", and thesis Chapter 2
(PDF pp. 16–33) IS that paper verbatim** — the 2001 data were inside a PDF already in the
scratchpad the whole time.
  Recovered (all 120 °C, 150 mM sugar, 3% Na-caseinate, 0.1 M phosphate, pH₀ 6.8): full
  glucose–casein and fructose–casein trajectories for glucose, fructose, formic acid, acetic
  acid, lysine residues and Amadori compound at 0/2/5/10/15/20/30/40 min; pH and total acid
  by titration AND by HPLC; A420 split total/protein/sugar plus protein-bound melanoidins;
  the isolated-Amadori control (Fig 2.7). **Plus the complete fitted kinetic model:** Scheme
  4.1's 9 ODEs, Table 4.1's 11 rate constants ± 95% CI at 120 °C, **Table 4.2's Arrhenius set
  at 90/100/110/120/130 °C with Ea per step**, and Schemes 4.3/4.4 with the pH correction and
  arginine. Tables 4.1/4.5 and all Scheme-4.4 ODEs verified digit-for-digit against 400-dpi
  renders.
  **Four independent validations** (none fitted): t=0 glucose reads 150.0 vs a 150 mM nominal
  charge; Fig 2.3's HPLC acids (different figure, different axis) match Fig 2.1/2.2's
  formic+acetic to within 0.12 mmol/L at all ten points; the thesis's own stated pH drops
  (0.3 and 0.4 units) match the digitized 6.73→6.45 and 6.71→6.30; mass balance closes at
  95% (101% including unidentified acids).
  **Source defects flagged as-is:** Table 4.3 (pH 5.9) is cited on thesis p.55 but is not
  printed anywhere; `d[M]/dt` uses `k14` where the system uses `k15` (thesis typo).
  Data in `<scratchpad>/brands_out_tables.md` (362 lines) with overlay PNGs. **[O]**

### Owner decisions this wave surfaces

1. **Retract the false "not extractable" note** in `brands_sugar_casein_120C_pH68.yml` and
   promote the trajectories as `data_quality: digitized_figure`. It is the reason a rich,
   already-downloaded dataset sat unused. The Chapter 4 rate constants are arguably the
   bigger prize: a published, fully-specified competing kinetic model at five temperatures
   for exactly this chemistry — a far stronger comparator than any single benchmark row. [O]
2. **Route `get_ph_multiplier` into the prediction path, or delete it.** Right now the repo
   has tested pH physics that no prediction ever sees, which is worse than having none —
   a reader inspecting `conditions.py` will conclude the model handles pH correctly. [O]
3. **Guard pH and aw recommendations at runtime.** 2/10 is worse than chance and the failure
   is systematic, not noisy. A user who moves pH on this model's advice is more likely to be
   misled than helped. [O]
4. **Trikusuma citation year → 2020**, and note the Processed-column provenance. [O]
5. **Hofmann: either get the body, or move the benchmark to mol %.** [O]

### What this score does and does not license (full text in the report)

It licenses **screening** claims of a narrow shape — "which arm gives more of X" with pH and
moisture held fixed — at roughly 84% on those axes. It licenses **nothing quantitative**: the
panel never scores a magnitude, and where magnitudes are visible they are often poor even
when the sign is right (SUG-13 gets glucose > fructose for furfural by 4.9x where the
measurement is 1.2x). It licenses **no pH or moisture recommendation at all**.
**62% is not a quality certificate** — chance is ~50%.

Two caveats on the panel itself: it is **scoped to the model's own footprint** (chemistry the
model cannot emit is recorded as unevaluable, not as failure, so this is not coverage of
Maillard chemistry); and several literature systems are food matrices mapped onto a
two-precursor aqueous model at assumed pH/aw, with the loose mappings flagged per claim.

**If this panel is re-run after a model change, report the strictly-independent number and
the per-category breakdown TOGETHER.** A headline that rose while pH stayed at 2/7 would be
hiding the finding.

---

## Wave S1 — the additive flux propagator + the unreachable matrix registry (2026-08-27)

Owner-approved. TWO STRUCTURAL FIXES in `src/recommend.py`, then one regeneration and
reconciliation. NO barrier, NO observability factor, NO projection constant was touched —
nothing was tuned to compensate for either fix, and the places where the numbers got worse
are reported as such. Fix 2 was landed and measured FIRST, then Fix 1, so the two effects
are separately attributed below.

### (a) THE ALGORITHM, BEFORE — stated so the change is reviewable

`Recommender.predict_from_steps` runs a Bellman-Ford-style relaxation over the enumerated
elementary steps. State per species: `tracking[canon] = (span, conc, depth, weight, unc)`.

  1. SPAN is a cumulative SERIES RESISTANCE, not a max: for a step of barrier `b` fired from
     reactants whose worst span is `s`, `exp(span/RT) = exp(s/RT) + exp(b/RT)` (log-sum-exp
     for numerical safety). The single largest barrier on a route dominates it.
  2. WEIGHT is the route's flux proxy:
     `path_weight = min_r_conc * co_reactant_factor * exp(-path_span/RT)`, times
     `_temporal_accessibility(tau*depth, t)` in FAST mode or an integrated Arrhenius
     propensity under a ramp.
  3. RELAXATION WAS WINNER-TAKES-ALL. `tracking[p]` was overwritten only `if path_span <
     current_span` (or, on an exact tie, if `path_weight > current_weight`). So at the
     fixpoint each product held the LOWEST-SPAN route and THAT ROUTE'S FLUX ALONE. Every
     other route to that product was discarded entirely.
  4. ALLOCATION (`_project_weighted_flux_to_ppb`, the Wave H budget machinery): the volatile
     targets that accumulated (`depth > 0`, finite span, in `target_lookup`, budget-relevant)
     get `activity = (weight/max_weight) * depth_activity * direct_sulfur_bonus`, then
     `mol_fraction = activity / total_activity` and
     `ppb = total_volatile_budget_molar * mol_fraction * MW * ppb_conversion_factor`. The
     budget itself comes from `_estimate_projection_budget` (conversion extent x limiting
     precursor x volatile yield fraction) and is INDEPENDENT of the propagator.

  CONSEQUENCE, which is Wave P's finding: because step 3 discarded every non-fastest route,
  step 4 saw one channel per product. Adding a second real route could not change a number
  unless it was FASTER than the incumbent.

### (a') THE ALGORITHM, AFTER

  Steps 1, 2 and 4 are UNCHANGED. Step 3 is unchanged for span/conc/depth/`best_paths`, so
  every span-driven behaviour downstream (span ranking, path traces, `depth_activity`) is
  bit-identical. What is new is ONE additional sweep, after the fixpoint:

      for step in steps:                       # from the CONVERGED tracking state
          ev = _evaluate_step(step)            # same arithmetic, extracted verbatim
          cid = _route_channel_id(ev.path)     # the route's FULL ORDERED STEP-SET
          for p in ev.products:
              channel_flux[p][cid] = max(channel_flux[p].get(cid, -inf), ev.path_weight)
      summed_channel_flux[p] = sum(channel_flux[p].values())

  and step 4 now reads `summed_channel_flux[p]` instead of `tracking[p][3]`. The sweep runs
  ONCE after convergence rather than accumulating during relaxation, which is what makes it
  impossible to add the same channel once per iteration. `_project_weighted_flux_to_ppb`
  keeps a `channel_flux_totals=None` path that reproduces the old winner-takes-all behaviour
  exactly; it is used by the mass-honesty test and by nothing that ships.
  The per-product channel breakdown is exposed as `debug_channel_flux` in the returned
  payload, so "the model has N routes to X" is now auditable against what was actually summed.

  THE REFACTOR IS PROVABLY INERT. `_evaluate_step` was extracted VERBATIM from the loop body,
  so the risk was a silent arithmetic change hiding inside a "pure" refactor. Measured, not
  assumed: with the additive layer disabled (`channel_flux_totals=None`) the live panel
  reproduces the pre-Wave-S1 predictions on 36 of 42 rows BIT-IDENTICALLY, and the 6 that
  differ are EXACTLY Fix 2's six rows. Against the Fix-2-only measurement it is 42 of 42
  identical. So every movement reported below is attributable to one of the two fixes and
  none of it to the extraction.

### (b) THE DEDUPE RULE — and the one the brief suggested, implemented, MEASURED and REJECTED

  SHIPPED RULE: two routes are the same channel iff they are the same ordered step-set.
  Within a channel the largest flux wins; across channels the fluxes SUM.

  REJECTED RULE: "two routes sharing their RATE-LIMITING STEP are one channel, take the max."
  It was implemented first, exactly as the brief specified, and then measured. TWO findings
  killed it, and the brief's own instruction was to verify before expecting the sum:

  1. THE RATE-LIMITING STEPS ARE NOT DISTINCT. Measured on cys_ribose_140C_Hofmann1998:
     BOTH MFT routes have the same highest-barrier step, `Amadori_Rearrangement` at 29.06
     kcal/mol. So does the FFT route. That step sits on the shared cysteine/ribose trunk that
     essentially every route in the network passes through:
         MFT winner path barriers  20.672 (Schiff) 29.060 (Amadori) 28.000 28.000 26.350
         FFT winner path barriers  20.672 (Schiff) 29.060 (Amadori) 21.000 28.000 23.300 26.800
     Under the rejected rule MFT keeps 242.38 ppb EXACTLY, FFT keeps 217.99, and the whole
     live panel moves 3 rows of 42 (resconi furfural x0.99992, Bolton MFT x1.0210, Cerny2008
     MFT x1.1464). That is winner-takes-all in all but name, and the [P] would not be closed.
  2. AND IT IS WRONG, not merely inert. X --(slow, R_c)--> Y, then Y --(fast)--> P by two
     branches and Y --(fast)--> Q by one. At steady state the trunk fixes the total flux and
     the branches PARTITION it by conductance: P's share is
     (1/R_1 + 1/R_2)/(1/R_1 + 1/R_2 + 1/R_3) = 2/3 for equal branches. This propagator's
     per-route weight is `pool * exp(-span/RT)` with `exp(span/RT) = sum_i exp(b_i/RT)`, so a
     dominant trunk collapses EVERY route's weight onto the same value ~ pool/exp(R_c/RT).
     SUMMING P's two branches then reproduces the 2/3; taking the max returns 1/2. A shared
     bottleneck is the reason the sum is correct, not a reason to suppress it — and the mass
     the sum appears to create is not created, because step 4 normalises against a fixed
     budget (verified in (d)).

  LIMITATIONS OF THE SHIPPED RULE, stated honestly and carried as [P]:
   * Every distinct producing step is treated as an independent branch. The model carries NO
     POOL DEPLETION, so two branches drawing on the same scarce intermediate are summed as if
     that intermediate were unlimited. Only the global budget cap stands between that and
     minted mass, and it bounds the TOTAL, not the split.
   * Because `best_paths` keeps one upstream route per reactant, a candidate's step-set is
     determined by its terminal step, so in the current enumeration the rule reduces to "one
     contribution per distinct producing step". A genuinely branched upstream is not enumerated.
   * Additivity is applied where flux is CONSUMED (at the product, by the allocation layer)
     and is NOT propagated. An intermediate with two routes still hands its own products a
     single span, so parallelism does not compound along a chain.

### (c) FIX 2 FIRST — the unreachable matrix registry (Wave O finding (f)), MEASURED ALONE

  CAUSE, located: `_apply_output_projection` called
  `describe_matrix_calibration(species.label or canon, ...)`. Species injected into
  `corrected_initial` by the lipid-oxidation path are materialised as
  `Species(species_name_lookup.get(canon, canon), canon)`, and when no enumerated step ever
  names them the label falls back to the CANONICAL SMILES. `describe_matrix_calibration`
  normalises "CCCCCC=O" to the literal string "ccccccc=o", matches no record and no class
  anchor, returns `calibration_observable_factor=None`, and the caller's `or 1.0` applied a
  factor of ONE, silently.

  FIX: `_registry_compound_name(canon, label, target_lookup)` at the lookup boundary — keep
  the species' own label when it IS a name, fall back to the compound-database name for that
  canonical SMILES when the label is merely the SMILES again. `target_lookup` is keyed by
  canonical SMILES and its `name` is the same spelling `MATRIX_BENCHMARK_BASE_MARKER_YIELDS`
  uses on the `matrix_only` lane, so the two lanes now agree on the key. Applied to BOTH
  branches of `_apply_output_projection`.

  MEASURED ALONE, against the Wave R tree, before Fix 1 existed:
    SCORED PANEL: 6 of 42 rows moved, all four internal snapshots, nothing else.
      pea_isolate_..._Internal2026 / Hexanal   0.171992 -> 0.742533   x4.317250
      pea_isolate_..._ProtocolPilot2026 / Hexanal      same          x4.317250
      soy_isolate_..._Internal2026 / Hexanal   0.178260 -> 1.700616   x9.540070
      soy_isolate_..._ProtocolPilot2026 / Hexanal      same          x9.540070
      soy_isolate_..._Internal2026 / Nonanal   0.039858 -> 0.042515   x1.066667
      soy_isolate_..._ProtocolPilot2026 / Nonanal      same          x1.066667
      pea Nonanal did NOT move: its ambient factor is 1.0.
      x4.317250 and x9.540070 are the Wave O refitted pea/soy ambient hexanal factors,
      finally applying. x1.066667 = 0.160/0.150, the soy-vs-pea ambient release ratio.
    EXTERNAL HOLD-OUT (n=200, seed 0): **ZERO of 8 points moved, artifact bit-identical.**
    THE PROOF WAVE O ASKED FOR, now a test: perturbing the reached record by 4.32x moves the
      soy Internal2026 Hexanal prediction by exactly 4.32x, and the factor is restored in the
      same test. NOTE WHICH RECORD IS REACHED, and it is Wave O's second observation:
      `determine_matrix_process_state(100 C, 45 min, aw 0.95)` returns `intermediate_matrix`,
      whose fallback order is (`intermediate_matrix`, `ambient_slurry`), so the record reached
      is the AMBIENT one — not the `heated_matrix` entry, and not the
      `aqueous_pre_extrusion_model` the benchmark file declares. That process-state mismatch
      is UNFIXED and remains open; this wave fixed the lookup, not the state.
    BEFORE THE SNAPSHOT REFRESH the four snapshots FAILED (pea max_ratio 1.0 -> 4.317,
      soy 1.0 -> 9.540, overall pass -> ranking-gap / scale-gap). They are model-generated
      reproducibility baselines, so the Wave M sequence refreshes them first; after the
      refresh the 4/4 internal-synthetic passes are restored. Recorded because the
      intermediate state is what a reviewer would see if they ran the panel before the refresh.

### (c') FIX 1 SECOND — the additive propagator, MEASURED ALONE on top of Fix 2

    SCORED PANEL: 23 of 42 rows moved.
      cys_ribose_140C_Hofmann1998 / 2-furfurylthiol        217.9999 -> 297.2755   x1.363650
      cys_ribose_140C_Hofmann1998 / 2-methyl-3-furanthiol  242.3782 -> 283.5889   x1.170027
      resconi_2023_pbma_beef_identity / furfural          3330.9660 -> 3148.6668  x0.945271
      thiamine_cys_xylose_145C_Cerny2008 / MFT              0.772987 -> 0.886115  x1.146351
      thiamine_cys_glucose_120C_Bolton1994 / MFT            0.017025 -> 0.017379  x1.020810
      the 4 internal snapshots, per compound:
          2-furfurylthiol                x1.499642
          2-methyl-3-furanthiol          x1.319008
          2,5-dimethylpyrazine           x0.974722
          bis(2-methyl-3-furyl) disulfide x0.974722
          furfural                       x0.974722
          Hexanal / Nonanal              unchanged (injected, not network-propagated)
      READ THE SIGNS TOGETHER: the sulfur channels RISE and the pyrazine/disulfide/furfural
      rows fall by EXACTLY the same factor to pay for them. That is the fixed budget
      redistributing, and it is the qualitative signature the mass-honesty check quantifies.
    EXTERNAL HOLD-OUT: **ZERO of 8 points moved, artifact bit-identical.**
    THE C2+C3 LANE NOW CONTRIBUTES. On Hofmann1998, with the lane monkey-patched off:
      pentodiulose lane alone 217.25 ppb, both lanes 283.59 ppb, i.e. +30.5%. Wave P's
      projected 242.38 + 71.02 = 313.39 was never obtainable: the lanes are not independent
      (shared rate-limiting step) AND the budget is fixed, so the single-lane figures cannot
      be added — the pentodiulose-alone figure itself falls 242.38 -> 217.25 once the
      competing FFT channels also become additive.
    CHANNEL CENSUS on Hofmann1998: 11 of 35 tracked species now carry >1 channel (was 8 under
      the rejected rule). MFT 2 channels, sum/max = 1.4006; FFT 2 channels, 1.6324.

### (d) MASS HONESTY — VERIFIED, NOT ASSERTED

  Method: intercept `_project_weighted_flux_to_ppb` and run the SAME converged state twice,
  once with `channel_flux_totals=None` (pre-Wave-S1 winner-takes-all) and once with the summed
  channels, then convert both allocations back to molar and compare against
  `projection_budget.total_volatile_budget_molar`.

    | system                        | budget (molar) | allocated/budget OLD | NEW  | summed ppb   |
    | ----------------------------- | -------------- | -------------------- | ---- | ------------ |
    | cys_ribose_140C_Hofmann1998   | 1.0442e-05     | 1.000000000000       | 1.000000000000 | 1152.26 -> 1158.83 (x1.0057) |
    | resconi_2023_pbma_beef        | 9.5114e-05     | 1.000000000000       | 1.000000000000 | 10848.47 -> 10849.06 (x1.00005) |
    | soy_..._Internal2026          | 7.6665e-09     | 1.000000000000       | 1.000000000000 | 0.7517 -> 0.7548 (x1.0042) |

  The molar total is the budget to 1 part in 10^12 before AND after. The summed PPB moves by
  <=0.6% and only because ppb is molecular-weight-weighted and the allocation shifted between
  species of different mass. Adding channels moves the SPLIT; it cannot mint mass.
  Pinned in `tests/scientific/test_wave_s1_additive_flux_2026_08.py::
  test_the_volatile_budget_still_caps_the_sum`.

### (e) THE FULL OLD -> NEW HEADLINE TABLE (Wave R tree -> Wave S1 tree)

  BENCHMARK PANEL (`benchmark_summary.json`), max_ratio / MALE:
    cys_ribose_140C_Hofmann1998          1.4110 / 0.0935 -> **1.4864 / 0.1267**  (Fix 1) WORSE
    resconi_2023_pbma_beef_identity      4.6573 / 0.6681 -> **4.4024 / 0.6437**  (Fix 1)
    thiamine_cys_glucose_120C_Bolton1994 763.5881/2.8829 -> **748.0216 / 2.8739** (Fix 1)
    thiamine_cys_xylose_145C_Cerny2008   3.1954 / 0.5045 -> **2.7874 / 0.4452**  (Fix 1)
    every other row bit-identical (the four snapshots were refreshed and recover exactly).
    THE HOFMANN CONTRACT (1.45x / 0.09 dex, UNTOUCHED) now fails on BOTH criteria where Wave
    P had it failing on one. MFT improved 1.4110x under -> 1.2060x under; FFT degraded
    1.0900x over -> 1.4864x over. Not clawed back: the lanes share their upstream trunk.
    status_counts scale-gap 8 / pass-no-ranking 2 / pass 4  ALL UNCHANGED.
    evidence-role split 6/4/4 UNCHANGED. predictive passes 0/6 UNCHANGED. fit-recovery 0/4
    UNCHANGED. internal-synthetic 4/4 UNCHANGED. Aggregate strict 4/14 UNCHANGED, lenient
    (presentation-layer, counts `pass-no-ranking`) 6/14 UNCHANGED. strict-ready 0/14 UNCHANGED.

  EXTERNAL HOLD-OUT (`external_validation_report.json`, n=200 seed 0) — ALL EIGHT POINTS AND
  EVERY SUMMARY FIELD **BIT-IDENTICAL**. The `cmp` over the artifact reports 0 leaf diffs.
      | point                          | measured | p50      | fold      | moved? |
      | bi_2020_raw_pea / hexanal      |   1260.0 |   1013.0 |    1.2437 | no |
      | bi_2020_roasted_pea / hexanal  |    324.0 | 801700.4 | 2474.3839 | no |
      | li_2026_hme / 1-hexanol        |    20.04 |  22394.3 | 1117.4799 | no |
      | li_2026_hme / 2-pentylfuran    |   5625.8 |  11038.6 |    1.9621 | no |
      | li_2026_hme / hexanal          |    605.6 |  56409.0 |   93.1471 | no |
      | li_2026_hme / nonanal          |    72.66 |   8596.4 |  118.3093 | no |
      | liu_2023_ppi / hexanal         |  11320.0 |   1013.0 |   11.1747 | no |
      | liu_2023_ppi / nonanal         |   0.8018 |    75.55 |   94.2231 | no |
    median |log10| **1.9717 dex UNCHANGED**; median fold **93.6837x UNCHANGED**;
    ci_coverage 3/8 (0.375) UNCHANGED; pre-widening 3/8 UNCHANGED;
    genuine_extrapolation 1/5 UNCHANGED; max_fold_error 2474.3839x UNCHANGED.
    THE REASON IS STRUCTURAL AND IS ITSELF A FINDING: all four hold-out bundles execute the
    `matrix_only` path. It passes compound NAMES to the registry (Fix 2 never applied) and it
    bypasses `predict_from_steps` entirely (Fix 1 never applied). The external hold-out
    exercises the lipid-oxidation and observability lane and says NOTHING about the Maillard
    network propagator — an eight-point insensitivity no previous wave had measured.

  MC PANEL (`prediction_uncertainty.json`, n=200 seed 0): every COUNT unchanged —
    benchmark_count 11, matched rows 35, aggregate coverage 29/35 (82.9%),
    honest_literature_coverage 1/3, not_evaluable 4, excluded_fitted_rows 2/2,
    excluded_..._would_have_been_hits 2. Only interval widths moved:
      external_literature median CI  0.8495 -> **0.9463** dex  (WIDER, coverage still 1/3)
      fitted_row                     2.2767 -> **2.2083** dex
      internal_synthetic             3.6929 -> **3.5612** dex
    THE EXTERNAL WIDENING IS HONEST AND IS NOT AN IMPROVEMENT: a compound reached by several
    routes now samples the barrier uncertainty of ALL of them, so the interval grew 1.25x and
    bought no coverage with it.

  VALIDATION OVERVIEW: benchmark_count 14, strict_ready 0, inside_1_5x 3, outside_1_5x 6,
    worst_quantitative_ratio 1203.6799x (CML) — ALL UNCHANGED. Only
    reference_worst_quantitative_ratio 3.1954 -> **2.7874** (Cerny 2008 MFT).
    family_validation_overview: SLR-09 carbohydrate pyrolysis mean |log10| 0.223 -> **0.215**;
    every other family row unchanged.

  MATRIX SIGMA (`derive_matrix_sigma_from_residuals.py`): artifact **BIT-IDENTICAL**.
    n=5, rms_ln_sigma 3.0166, bias_fold 3.3149, centered_sd 3.0951, 90% CI [2.0274, 6.3038],
    shipped 2.86 INSIDE and NOT moved. Structural, per Wave O's own note: the uncalibrated
    tier multiplies oxidation load by a base marker yield and reads NEITHER an observability
    factor NOR the propagator. Neither fix can reach it, in either direction.

  PROJECTION REFIT (`refit_projection_constants.py`, still NOT applied): fitted tau UNCHANGED
    at 10000 min; objective at the fit 0.8817 -> **0.8824** dex, at the shipped tau 12589
    0.8879 -> **0.8863** dex. The budget scale the fit wants is still 1.26x.

  SNAPSHOTS: `refresh_internal_reproducibility_snapshots.py` GENERATOR_TAG v8 -> **v9** with a
    new leading Wave S1 block that explicitly marks the Wave P "EXACTLY ZERO" paragraph and
    the Wave O "THIS SNAPSHOT DID NOT MOVE AT ALL / registry UNREACHABLE" paragraph as
    SUPERSEDED, rather than deleting them. NOTE: Wave P's ledger entry claims it bumped
    "v8 -> v9", but the shipped tag in the four snapshot files was still v8 — Wave P edited
    the NOTES text without the tag. The tag is now genuinely v9. [P] minor, recorded.

  PENTOSE >> HEXOSE: 6.15x -> **7.78x** (ribose 686.83 -> 824.72, glucose 111.65 -> 105.95).
    NOT IMPROVED SUGAR DISCRIMINATION, and the guard test says so in its own failure message.
    Structural share re-measured in-process (thiol_addition_pentodiulose set equal to
    thiol_addition_hexose): **3.1368x**, up from 2.3118x. History of the split: 1.13x of
    3.39x (Wave N), 2.31x of 6.15x (Wave P), 3.14x of 7.78x now. The structural share HAS
    grown, for a defensible reason — the propagator stopped discarding the pentose limb's
    parallel routes — but ~2.5x of the claim is still the gap between a FITTED barrier and an
    UNCONSTRAINED LEGACY FIT. Both halves are in the README and AUDIT.

  MECHANISTIC-PRIORITY BENCHMARKS: 2 -> **0**, and this is a LOSS OF SIGNAL, not a win.
    Fix 2 gives the internal snapshots' Hexanal/Nonanal rows a real compound-specific record,
    so they stop reading `calibration_evidence_strength = "heuristic"` and
    `_matrix_closure_action` stops routing them to `mechanistic_blocker`. What they read
    instead is `process_state_mismatch` — the registry's honest label for a record reached
    only through the `intermediate_matrix -> ambient_slurry` fallback — and the rule set has
    NO BRANCH FOR THAT LABEL, so it is scored exactly like a genuine class-anchored transfer.
    The governing decision did NOT change (`hold_observable_first`, 0 approved offline jobs,
    same three blockers), so no compute was unlocked; what was lost is a warning. Pinned at 0
    with the cause in three test files; the rule-set repair is [P] and was deliberately NOT
    made in the same pass as the fix that exposed it.

### (f) TESTS RE-PINNED (each with a dated causal comment; none relaxed)

  Entry state after the two fixes + regeneration: 8 failures, all in tests/scientific.
  tests/unit + tests/integration + tests/scripts were green throughout (1056 passed,
  1 skipped, 0 failed, measured before any re-pin).

  1. `test_wave_p_chemistry_2026_08.py::test_the_second_mft_channel_contributes_exactly_zero_
     to_the_prediction` — WAVE P'S OWN FLAG, and its failure message said this failure would
     be the notification. RENAMED to `..._now_contributes_after_the_wave_s1_propagator_fix`
     and INVERTED: pentodiulose-alone 217.25, both lanes 283.59, ratio pinned at 1.3054x. The
     docstring carries Wave P's original text verbatim plus why the 313.39 projection never
     applied.
  2. `test_wave_p_chemistry_2026_08.py::test_hofmann1998_after_the_refit` — MFT 242.38 ->
     283.59, FFT 217.99 -> 297.28, with the full "half of this got worse" note and the
     statement that no barrier was refitted.
  3. `test_free_aa_quantitative_regression.py` — MFT band (1.15, 1.411, 1.75) ->
     (1.00, 1.206, 1.50); FFT (1.00, 1.090, 1.43) -> (1.00, 1.486, 1.95). Upper bounds carry
     the Wave P pins' RELATIVE spans (x1.240 and x1.3119); the MFT lower bound is CLIPPED at
     1.00 because a symmetric fold error cannot go below it, which makes that side STRICTER,
     not looser. The header block now records that the contract fails on BOTH criteria again.
  4. `test_honest_headline_guards.py::test_honest_external_literature_coverage_is_1_of_3_...`
     — median CI width 0.8495 -> 0.9463 dex, with the companion widths and an explicit note
     that a widening which buys no coverage is not an improvement. All six COUNT assertions
     in the same test are untouched and still pass, which is what identifies this as an
     interval-width change.
  5. `test_honest_headline_guards.py::test_pentose_hexose_mft_ordering_is_6_15x_...` —
     RENAMED to `..._is_7_78x_not_the_retired_6_15x_3_39x_8_98x_or_15_8x`; ribose 686.8 ->
     824.7, glucose 111.6 -> 106.0, ratio 6.15 -> 7.78, structural share 2.31 -> 3.14, and the
     README/AUDIT doc-token assertion moved 6.15 -> 7.78 (both docs updated in the same edit).
  6. `test_matrix_observable_closure_audit.py` — Hexanal closure_action `mechanistic_blocker`
     -> `class_level_transfer_acceptable`; watchlist -> empty. TIGHTENED while re-pinning: it
     now ALSO asserts `calibration_evidence_strength == "process_state_mismatch"` (so the
     movement's actual cause is pinned, not just its effect) and pins the renderer's
     "Mechanistic watchlist count: 0" line rather than the mere presence of a heading.
  7. `test_offline_refinement_governance.py` — `mechanistic_priority_benchmark_count >= 1`
     -> `== 0`. TIGHTENED: the assertion was a FLOOR, which could not have detected movement
     in either direction; it is now two-sided, and the priority list is pinned empty with a
     failure message naming the SMILES-vs-name regression as the thing to look for.
  8. `test_refinement_governance.py` — same count 2 -> 0, same treatment. Both tests keep
     their `governing_status == "hold_observable_first"` and
     `approved_offline_job_count == 0` assertions untouched, which is the evidence that no
     compute was unlocked by the change.
  9. NEW `tests/scientific/test_wave_s1_additive_flux_2026_08.py` (9 tests) — two synthetic
     independent routes to one product ADD exactly (and carry measurably different flux, so
     max and sum cannot agree); a duplicated route contributes ONCE; channel identity is the
     step-set and not the rate-limiting step; both Hofmann MFT routes are pinned as sharing
     `Amadori_Rearrangement`; MASS HONESTY at a fixed budget; the Hofmann pair after the fix;
     the Wave O 4.32x perturbation now moving the prediction; `_registry_compound_name`'s four
     cases; and the exact set of rows Fix 2 moved.

  NOT RE-PINNED, DELIBERATELY: `cys_ribose_140C_Hofmann1998`'s own contract (1.45x / 0.09 dex,
  now failing on both criteria), the 3.0x pentose/hexose floor in
  `test_pentose_hexose_sulfur_ordering.py`, the Pratap-Singh ratio tolerance, the Trikusuma
  observability factors, `max_fold_error`, the pre-widening 1/5, the 0/6 predictive headline,
  the 4/14 and 6/14 aggregates, strict-ready 0/14.

### (g) DOC SYNC

  README.md: the sulfur trust paragraph now reads 1.21x under on MFT AND names the 1.49x FFT
    over-prediction as worse than before; pentose block 6.15x -> 7.78x with the structural
    share 2.31x -> 3.14x and the honest "the structural share has grown but a third of the
    claim is still a barrier gap"; the Wave P refit block's closing sentence now says its two
    numbers have since moved and got worse; a NEW "Wave S1" block with the old->new table for
    MFT/FFT/max_ratio/MALE, the rejected-dedupe finding, the correction of Wave P's 313.39
    arithmetic, and the mass-honesty verification.
  AUDIT.md: the "what survives is ordering" headline 6.15x -> 7.78x with the full split
    history; 14 new rows in the "What Round 3 cost" table, every one labelled *(Wave S1)*,
    including the three that got worse and the one that is a loss of signal; a NEW section
    "Wave S1 — the additive flux propagator, and the registry nobody could reach" with the
    per-fix attribution table, the rejected rule and its two-part refutation, the mass-honesty
    numbers, and the hold-out insensitivity finding; and the Wave P "EXACTLY ZERO" paragraph
    now carries a "-> CLOSED by Wave S1" marker instead of standing as current.
  docs/reference/VALIDATION_CONTRACT.md §3E: a new dated paragraph recording that Wave S1
    moved 26 of 42 in-panel rows and ZERO of the 8 hold-out points, with the structural reason
    (`matrix_only` bypasses both fixes) and the conclusion that the hold-out's invariance is
    evidence about the hold-out's coverage, not about the model.

### (h) [P] CARRIED FORWARD

  1. The shipped dedupe rule treats every distinct producing step as an independent branch,
     and the model has NO POOL DEPLETION. Two branches drawing on the same scarce intermediate
     are summed as if it were unlimited; only the global budget cap bounds the result, and it
     bounds the total rather than the split.
  2. Additivity is not PROPAGATED: an intermediate with two routes still hands its products a
     single span, so parallelism does not compound along a chain. A full fix means carrying
     summed flux through the relaxation, which is a larger rewrite.
  3. `_matrix_closure_action` has no branch for `process_state_mismatch`, so a calibration
     factor reached only through a process-state fallback now scores as an acceptable
     class-level transfer. This is what silently emptied the mechanistic-refinement watchlist.
  4. The runtime recomputes the internal snapshots' process_state as `intermediate_matrix`
     from 100 C / 45 min / aw 0.95, not the `aqueous_pre_extrusion_model` those benchmark
     files declare. Wave O flagged it; still open. It is why the AMBIENT record, not the
     heated one, is the record Fix 2 made reachable.
  5. The external hold-out cannot see the Maillard network at all — all four bundles are
     `matrix_only`. Until a hold-out bundle exercises `predict_from_steps`, no external
     evidence bears on the propagator, the barriers, or the network topology.
  6. Wave P's ledger claims a snapshot NOTES bump to v9 that never reached GENERATOR_TAG
     (the files shipped v8). Minor, and now corrected.
  7. Everything Wave P carried forward that this wave did not touch: the `src/conditions.py`
     family-name substring coupling; the C2+C3 lane's missing moisture dependence; the oleate
     pool on the linoleate rate constant; 2-mercapto-3-pentanone's missing odour threshold;
     the two disagreeing panel pass-counters (4/14 strict vs 6/14 presentation-layer); the
     Trikusuma pre-Wave-P MALE inconsistency; the isotope dossier's Hearn & Smith
     mis-attribution; the `[HH]` gate still on `_thiol_addition` (FFT); and the missing hexose
     retro-aldol route to glycolaldehyde.

### (i) GATES + SUITE

  GATES, all three re-run on the final tree, after every code, artifact, test and doc edit:
    citation_gate   PASS — 82 files, 964 DOI-bearing fields, 317 unique DOIs (was 81/951/316
                    at Wave P; the +1/+13/+1 is Wave R's, not this wave's — Wave S1 introduced
                    NO new citation). WAIVERS and TEXT_SURFACE_WAIVERS are both still the
                    empty tuple.
    holdout_guard   PASS — 3/3 invariants, re-run after the regeneration. Worth stating
                    plainly given (e): the hold-out was not merely excluded from every fit,
                    it was not even REACHABLE by either fix.
    fit_target_gate PASS — both checks. The fit-target set is unchanged; this wave declared
                    no new fit and moved no fitted constant.

  FIRST FULL SUITE (`tests/unit tests/scientific tests/integration tests/scripts`, documented
  conda path): **1274 passed, 1 skipped, 2 xfailed, 0 FAILED** in 877.27 s, exit code 0.
  Arithmetic: 1274 + 1 + 2 = **1277** = Wave P's certified 1268 + the 9 new tests in
  `tests/scientific/test_wave_s1_additive_flux_2026_08.py`. The 1 skip and both xfails are the
  declared, strict-marked ones from Wave J2 (`xfail_strict = true`, so neither can silently
  start passing).
  THAT RUN IS SUPERSEDED AND HERE IS WHY, because the honest thing is to say it rather than
  quote it: two edits landed while it was in flight — a one-key addition to the matrix-only
  payload in `src/benchmark_validation.py` (`debug_channel_flux: {}`, so both execution paths
  return the same key set) and the README/AUDIT prose sync. Python imports at collection, so
  neither was in effect for that run. The CERTIFYING run below was executed on the final tree
  with nothing further edited afterwards.

  CERTIFYING FULL SUITE, on the final tree (29 modified tracked files + the one new test file),
  after all three gates were re-run green on that same tree:
  **1274 passed, 1 skipped, 2 xfailed, 0 FAILED** in 871.98 s, exit code 0, zero FAILED/ERROR
  lines. It reproduces the first run's counts EXACTLY, which is the evidence that the two
  in-flight edits were inert. NOT COMMITTED, NOT STASHED — handed to the orchestrator as
  instructed.

  INTERMEDIATE RUNS, kept because they isolate the blast radius:
    - tests/unit + tests/integration + tests/scripts, after both fixes and the regeneration
      but before any re-pin: **1056 passed, 1 skipped, 0 failed** (417.22 s). Every failure
      this wave produced was in tests/scientific.
    - tests/scientific at that same point: **8 failed, 210 passed, 2 xfailed** (468.96 s).
      All 8 are the re-pins itemised in (f); none was a code defect.

  ENVIRONMENT NOTE. Free RAM fell to ~93 MB during the long runs and pytest advanced a few
  tests per minute through the tests/scientific governance region. Disk had 31 GiB free
  throughout, so this is Wave P's memory-pressure condition and NOT its ENOSPC condition. It
  slows the suite; it did not produce a single spurious failure in any run this wave.

---

## Wave S1b — the pH / water-activity ROUTING repair (2026-08-27)

Owner-approved. THREE ROUTING FIXES in `src/conditions.py`, then one regeneration and
reconciliation. **NO correction curve was reshaped, NO constant was refitted, NO barrier was
touched, and the directional panel was NOT iterated against** — the wiring was chosen from the
chemistry, measured once, and reported whichever way it came out. Four benchmark rows got
worse and are reported as such. Two data-curation notes were added; no measured value, DOI or
contract threshold was edited.

### (a) THE BASELINE, corrected before anything else

Wave S2 scored 18/29 on `HEAD = c1a12d2`. Wave S1 landed after that. **Re-running the identical
panel on `263bae8`, before any edit, gives 19/29** — Wave S1 had already moved SUG-04 (a
*fit_adjacent* row) from agree to `flat` and moved one independent row into agreement, so the
fit-adjacent bucket was 8/9 not 9/9 before this wave began. **19/29 is the honest baseline.**
pH and aw were untouched by Wave S1 and stood at exactly Wave S2's 2/7 and 0/3.

### (b) EACH DEFECT CONFIRMED BEFORE IT WAS TOUCHED

1. **`get_ph_multiplier` dead on the prediction path — CONFIRMED.** `grep -rn` finds callers
   only in `kinetics.py`, `pathway_ranker.py`, `cantera_export.py`; none is reachable from
   `evaluate_benchmark_payload`, which enters `conditions.get_rate_constant()` at
   `benchmark_validation.py:662`. That call site is the honest wiring point: it is where the
   per-family rate constant is computed and immediately converted to `bar_eff` for the
   propagator.
2. **Pyrazine branch unreachable — CONFIRMED, and quantified.** Enumerated every family the
   engine emits over all of `data/benchmarks/`: **29 distinct families, none containing
   "pyrazine"**. The pyrazine step is `Aminoketone_Condensation` (`reaction_templates.py:1000`).
3. **`_water_activity_correction` reaches almost nothing — CONFIRMED.** 3 of 29
   (`Amadori_Rearrangement`, `Strecker_Degradation`, `Lipid_Strecker_Synergy`), none on the
   furan/HMF track; its dehydration branch keyed on `"furfural"`, which matches no emitted
   family either, so that branch was dead too.

### (c) THE DESIGN QUESTION THE BRIEF ASKED — are the two pH terms the same physics?

**No. They are different terms, and that is why the sets are disjoint by construction.**

* `_ionization_correction` = **reagent availability**. Henderson-Hasselbalch fraction of the
  nucleophile present as free base. Monotone increasing in pH. Belongs to families where a
  nitrogen lone pair attacks.
* `get_ph_multiplier` = **route selection**. Which way the Amadori compound enolises:
  1,2-enolisation is acid-catalysed and opens 3-deoxyosone → furfural/HMF; 2,3-enolisation is
  base-catalysed and opens 1-deoxyosone → reductone → pyrazine/furanone. A branch-point
  partition, not a reagent count.

  `_ALPHA_AMINO_NUCLEOPHILE_FAMILIES` (pKa 8.0) — `Schiff_Base_Formation`, `Lipid_Schiff_Base`,
  `Amadori_Rearrangement`, `Strecker_Degradation`, `Lipid_Strecker_Synergy`. This reproduces
  EXACTLY what the old "amadori"/"strecker"/"schiff" substrings reached among emitted
  families: a de-substring-ing, not a widening.
  `_AMINOKETONE_NUCLEOPHILE_FAMILIES` (pKa 6.5, unchanged) — `Aminoketone_Condensation`.
  `_ENOLISATION_ROUTE_PH_FAMILIES` — `Enolisation_1_2`, `Enolisation_2_3`,
  `Enolisation_2_3_Amadori`. Disjointness is asserted at import AND behaviourally over all 29
  emitted families in the new test module.

`get_ph_multiplier` is deliberately NOT wired to the families its own substrings would sweep
up. `"thiol"`/`"thio"`/`"furan"`/`"cysteine"` reach the whole downstream sulfur and furan
track, and re-applying a 4.9x acid boost at every downstream addition, cyclodehydration and
oxidation would compound one physical effect five or six times along a route. And
`"condensation"` matches `Aminoketone_Condensation`: an ungated call would hand the PYRAZINE
step the ACID-peaked Schiff Gaussian, i.e. defect 2 again in the opposite direction. That is
what `_enolisation_route_ph_correction` exists to prevent, and it is pinned.

`Enolisation_Intermediate` is EXCLUDED: it is the step BEFORE the branch point, common to both
arms, so a factor there would scale both arms identically and express no selection.
`Fructofuranosyl_Dehydration` is EXCLUDED: it bypasses the Amadori branch point entirely. **That
second exclusion costs a panel row and is carried as an owner decision — see (g).**

### (d) THE aw MEMBERSHIP RULE — measured, not asserted

Criterion: **net water produced per step, counted directly off the enumerated steps** over all
of `data/benchmarks/`. Water-releasing → the (unchanged) `1.3 − aw` dehydration curve;
water-consuming → mass action in `aw`, floored the same way; net-zero → nothing.
`Amadori_Rearrangement` / `Strecker_Degradation` / `Lipid_Strecker_Synergy` keep the EMPIRICAL
Labuza peaked curve they already had, untouched — that curve is not a mass-action term, it is
the measured overall-browning response.

Every family is stoichiometrically UNIFORM except `Additive_Thermal_Degradation`
(+2 / 0 / −1 / −2), which is therefore in neither set: a single family-level factor cannot
honestly represent it. In the glycine/glucose system Wave S2 measured, the correction now
reaches **5 of 7** emitted families, up from 2.

LIMITATIONS, stated rather than hidden and both carried as [P]:
 * The inhibition applies ONCE per step regardless of whether the step sheds one water or
   three — the shipped `1.3 − aw` shape carries no stoichiometric exponent and adding one
   would be retuning the curve. A three-water dehydration is under-penalised relative to a
   one-water one.
 * For families that lump a redox into the same step (`Deoxyosone_Reduction`,
   `Mercaptoketone_Formation`, `Thiol_Addition_H2`, `Furanone_Reductive_Opening`,
   `Furan_Ring_Aromatisation`, and `Aminoketone_Condensation`'s H2) the released water is
   partly a redox byproduct, so the Le Chatelier argument is weaker for them. They are included
   anyway, because excluding them would need exactly the per-family judgement call the
   measurement replaced — but a reader should know the argument is weaker there.

### (e) TWO BUGS THIS WAVE INTRODUCED AND FOUND, both reported

1. **The Wave H substring guard short-circuited the furan track.**
   `_releases_rather_than_attacks_with_the_amine` matches `"enolisation"`, and it sat in front
   of the new water-activity sets, so `Enolisation_1_2` — the 3-deoxyosone → furfural/HMF
   dehydration, THE furan track, the one family finding 3 most needed to reach — returned 1.0
   before any set was consulted. Found by inspecting the factor table at the snapshot
   conditions, NOT by the score. The guard was removed from that function; the explicit sets
   subsume it strictly better (`Enolisation_2_3_Amadori` is net-zero, is in no set, still
   returns 1.0, and `tests/unit/test_wave_h_2026_08.py` still passes untouched). Every number
   in this entry is post-fix; the intermediate measurement (20/29 → wrong, headline 18/29) is
   superseded and is recorded in the report addendum so the sequence is auditable.
2. **`Furanone_Reductive_Opening` (net +1) was omitted from the set by transcription**, against
   this wave's own stated criterion. Caught by the new test that RE-DERIVES membership from the
   enumerated steps rather than trusting the source literal. It now checks all 29 emitted
   families and reports zero mismatches. Correcting it changed **no panel outcome and no
   benchmark row** — only DMHF's aw response.

### (f) THE PANEL, PRE → POST (measured once on the final wiring)

| Bucket | Wave S2 (`c1a12d2`) | Pre-fix (`263bae8`) | **Post-fix** |
|---|---|---|---|
| Strictly independent (headline) | 18/29 | **19/29 (66%)** | **20/29 (69%)** |
| + system-overlap | 6/6 | 6/6 | 6/6 |
| Screening total | 24/35 | 25/35 | **26/35** |
| Fit-adjacent (excluded) | 9/9 | 8/9 | 8/9 |

| Category | pre | post |
|---|---|---|
| **ph** | **2/7** | **4/7** |
| **moisture_aw** | **0/3** | **0/3** |
| **pH + aw** | **2/10 (20%)** | **4/10 (40%)** |
| sugar_identity | 8/8 | **7/8** |
| Excluding pH and aw | 17/19 (89%) | **16/19 (84%)** |
| additive_cysteine / temperature / time / lipid_lane / matrix_identity | 4/4 · 6/8 · 2/2 · 2/2 · 1/1 | unchanged |

PER-CLAIM FLIPS — every row whose outcome or shape moved:

| Claim | Observable | Expected | pre | post | |
|---|---|---|---|---|---|
| **PH-04** | 2,5-DMP, pH 4.5→5.6→6.5 | increasing | `decreasing` 64.5/38.8/24.4 | **`increasing`** 0.026/0.131/2.30 | **GAINED** |
| **PH-06** | 2,5-DMP, pH 4→7→9 | increasing | `decreasing` 99.5/16.8/14.9 | **`increasing`** 1.80/7.95/20.7 | **GAINED** |
| **SUG-12** | HMF, fructose vs glucose | A>B | `A>B` 2433 / 1524 | `A<B` 897 / 1617 | **LOST** |
| PH-07 | furfural, pH 4→7→9 | flat | `increasing` 752/902/908 | `peak` 793/880/792 | miss, closer |
| PH-05 | furfural, pH 4→5.5→7 | decreasing | `increasing` 305/681/894 | `increasing` 372/852/819 | miss |
| PH-03 | FFT, pH 4→5.5→7 | decreasing | `peak` 0/0.055/0.045 | `peak` 0/0.018/0.005 | miss |
| AW-01 | HMF vs aw | decreasing | `flat` 580/597/593 | `flat` 613/635/634 | miss |
| AW-02 | acrylamide vs aw | peak | `trough` 8.04/3.90/4.36 | `decreasing` 3.37/0.69/0.54 | miss |
| AW-03 | HMF vs aw | peak | `flat` 570/598/594 | `increasing` 600/636/634 | miss |
| SUG-04 *(fit-adjacent)* | MFT vs FFT | A>B | `flat` 284/297 | `A<B` 155/267 | already missing, now inverted |

**The two gained rows are exactly the two claims defect 2 was diagnosed against** — the two
independent direct measurements of dimethylpyrazine vs pH.

### (g) THE ONE LOST ROW, and the owner decision behind it [P]

SUG-12 is HMF from fructose vs glucose at pH 6.0. Glucose reaches HMF through
`Enolisation_1_2`, which now takes a **3.0x** acid boost at pH 6.0. Fructose reaches HMF
through `Fructofuranosyl_Dehydration`, excluded from that set, so it gets none, and glucose
overtakes.

**The argument runs both ways.** *For inclusion:* fructofuranosyl dehydration to HMF is the
textbook acid-catalysed hexose dehydration, and two routes to the SAME product now get
opposite pH treatment — an internal inconsistency. *For exclusion (shipped):*
`get_ph_multiplier`'s branch is documented as the 1,2- vs 2,3-ENOLISATION partition, and this
family is not part of that partition. **The exclusion was written down before the panel was
re-scored and was deliberately not revisited afterwards** — re-wiring it having seen that it
costs a row is precisely the optimisation-against-the-panel the brief forbade. **Owner decides
on the chemistry.**

### (h) WHY aw IS STILL 0/3 — structural, not routing

AW-01 and AW-03 both score **HMF** in a glycine/glucose system. The correction now reaches
that chemistry and visibly moves compounds — DMHF **115.4 → 45.5 → 49.6 ppb** over aw
0.25 → 0.65 → 0.95, a 2.5x separation where there was none; DMP 2.15 → 0.29 → 0.26. But **HMF
cannot respond**: HMF and furfural are both products of `Enolisation_1_2`, the SAME family, so
they always carry the same aw and pH factor, and together they are **90–96%** of that system's
volatile budget. Their shares are pinned against each other and the projection budget itself is
aw-independent. **No family-level correction can move HMF against aw in a two-precursor
system.** That is a statement about the allocation layer. [P]

AW-02 needs the acrylamide DESTRUCTION term Wave S2 identified as the model's most
consequential single miss; this wave did not add one.

### (i) THE COST — four benchmark rows, none clawed back

`status_counts` UNCHANGED (scale-gap 8 / pass-no-ranking 2 / pass 4), strict-ready UNCHANGED
0/14, panel size 14 UNCHANGED. The four internal snapshots recover exactly after the documented
refresh. Four real rows moved, all four the wrong way:

| Benchmark | max_ratio | MALE |
|---|---|---|
| `cys_ribose_140C_Hofmann1998` | 1.4864 → **2.2086** | 0.1267 → **0.2352** |
| `thiamine_cys_glucose_120C_Bolton1994` | 748.02 → **6730.85** | 2.8739 → **3.8281** |
| `thiamine_cys_xylose_145C_Cerny2008` | 2.7874 → **23.4061** | 0.4452 → **1.3693** |
| `resconi_2023_pbma_beef_identity_benchmark` | 4.4024 → **5.4570** | 0.6437 → **0.7370** |

ATTRIBUTION, measured one fix at a time by emptying the other two family sets at runtime (no
source edit, no constant moved):

| System / observable | none | fix 1 only | fix 2 only | fix 3 only | all three |
|---|---|---|---|---|---|
| Hofmann MFT (meas 342) | 283.6 | 219.7 | 296.5 | 286.2 | **154.8** |
| Hofmann FFT (meas 200) | 297.3 | 317.6 | 310.9 | 245.2 | **267.5** |
| Bolton MFT (meas 13.0) | 0.01738 | 0.01659 | 0.01744 | **0.002055** | **0.001931** |
| Resconi furfural | 3149 | 3721 | 3265 | 2922 | **3903** |

* **Bolton and Cerny are carried almost entirely by fix 3 (aw).** At aw 0.98 every
  water-shedding step takes `1.3 − 0.98 = 0.32`, i.e. **+0.89 kcal/mol** of effective barrier;
  the thiamine → MFT lane's terminal `Furan_Ring_Aromatisation` takes it while the
  `Additive_Thermal_Degradation` steps above it (non-uniform, excluded) do not. MFT is a tiny
  share of a large budget there, so its competitors absorb what it loses and the ratio amplifies.
* **Hofmann is carried by fix 1 (route pH), with interaction.** At pH 5.0 `Enolisation_1_2`
  gets 4.5x and `Enolisation_2_3_Amadori` ≈1.0, so the FFT/furfural arm gains on the MFT arm.
  Wave S1b's first draft called that "the textbook enolisation pH physics disagreeing with
  the benchmark — one of the two is wrong". **THAT WAS WRONG, and Wave S2b (same day, section
  below) is why.** There is no measurement on the other side: 342 / 200 ppb were derived
  inside this repository from `data/benchmarks/maillard_validation_benchmarks.md` §1.3, an
  abstract-reconstructed range table committed in the SAME commit as the benchmark file, and
  both are interior points of two OVERLAPPING invented bands (MFT 228-571 ppb, FFT
  114-342 ppb). The MFT > FFT ordering is midpoint selection; the 1.45x / 0.09 dex contract is
  ~1.7x tighter than its own source band. The mechanism and the degradation are both real; the
  yardstick is not. Corrected in README, AUDIT, the report addendum and two test comments, and
  the benchmark file's `wave_s2_followup` note now carries Wave S2b's finding plus its
  off-by-one correction (Cerny 2015's 145 C / 20 min sentence cites Hofmann & Schieberle
  **1995**, not 1998). Wave S2b's structural recommendations — retire the contract, demote out
  of PRIMARY, mark `no_verifiable_source` — are OWNER decisions and Wave S1b did not take them,
  and did not relax the contract either. See (l) and `## Wave S2b`.

OTHER REGENERATED HEADLINES:
  MC PANEL (n=200, seed 0): benchmark_count 11, matched rows 35, excluded_fitted_rows 2/2 —
    UNCHANGED. Coverage 29/35 → **28/35**. **`honest_literature_coverage` 1/3 → 0/3**: the one
    literature hit was Resconi furfural, whose `inside_ci` flipped True → False — its p50 rose
    2504.50 → 3462.43 ppb while its 90% CI NARROWED 1.4048 → 0.7460 dex. A narrower interval
    that now misses every literature row is the worst of both, and that is how it is reported.
    Companion widths: fitted_row 2.2083 → 2.9573, internal_synthetic 3.5612 → 4.1780.
  EXTERNAL HOLD-OUT: **all eight points and every summary field BIT-IDENTICAL** — median
    93.68x, coverage 3/8, max_fold_error 2474x, pre-widening 1/5, genuine extrapolation 1/5.
    Same structural reason as Wave S1: all four bundles run `matrix_only`, which never reaches
    `get_rate_constant`. **Two consecutive waves of prediction-path changes have now been
    invisible to the hold-out.** [P]
  VALIDATION OVERVIEW: inside_1.5x 3 → **2**, outside 6 → **7**; worst experimental point
    1203.680x (CML) → **6730.854x (Bolton MFT)**; worst reference-only 2.787x → **23.406x**;
    median authoritative compound max-ratio 108.218 → **112.227**. benchmark_count 14 and
    strict_ready 0 UNCHANGED. family_validation_overview: SLR-01 mean |log10| 0.179 → 0.283,
    SLR-09 0.215 → 0.246; every other family row unchanged.
  MATRIX SIGMA: artifact **BIT-IDENTICAL** — the derivation reads neither the propagator nor
    these corrections.
  PENTOSE ≫ HEXOSE: 7.78x → **18.27x**, structural share 3.14x → **4.27x**. **NOT improved
    sugar discrimination** and the guard test says so: BOTH numbers fell (ribose 824.7 → 374.0,
    glucose 106.0 → 20.5) because the hexose limb's `Thiol_Addition_Hexose_Legacy_Shortcut`
    sheds THREE waters and is penalised harder at aw 0.98 than the pentose limb's two steps.
    The ratio rose because the denominator fell faster. History of the split: 1.13x of 3.39x
    (N), 2.31x of 6.15x (P), 3.14x of 7.78x (S1), 4.27x of 18.27x now.
  SNAPSHOTS: `refresh_internal_reproducibility_snapshots.py` GENERATOR_TAG v9 → **v10** with a
    leading Wave S1b block. Movement, pea and soy identically — DMP x0.0022206 (penalised
    twice: 89% protonated at pH 5.6 AND sheds two waters), MFT x0.337603, disulfide x0.496265,
    FFT x0.504176, furfural x1.074000, Hexanal/Nonanal x1 (injected, never propagated).
    **THE RANKING CONTRACT CHANGED and it is a real movement, not a refresh artifact:** FFT now
    outranks MFT and the disulfide now outranks DMP.

### (j) THE DEAD-KEY CENSUS — seven, reported not silently deleted

Measured against the 41 family names the engine can emit:

| Key | Where | |
|---|---|---|
| `"pyrazine"` | `_ionization_correction`, `_water_activity_correction`, `get_ph_multiplier` | **DEAD** — fixed in the first two |
| `"furfural"` | `_water_activity_correction` | **DEAD** — this was the whole dehydration branch |
| `"heyns"` | `get_ph_multiplier` | **DEAD** |
| `"nitrogen_heterocycle"` | `get_ph_multiplier` | **DEAD** |
| `"oxygen_heterocycle"` | `get_ph_multiplier` | **DEAD** |
| `"1,2"` | `get_ph_multiplier` | **DEAD** — families spell it `1_2` |
| `"2,3"` | `get_ph_multiplier` | **DEAD** — families spell it `2_3` |

And one that is LIVE and was worse than dead: `"condensation"` in the Schiff branch matches
`Aminoketone_Condensation`. All seven are pinned as still-dead by
`test_the_substring_keys_that_matched_nothing_stay_documented_as_dead`, so a future rename into
one of them fails a test instead of moving a prediction.

### (k) TESTS RE-PINNED — each with a dated causal comment; none relaxed

Entry state: `tests/unit` + `tests/integration` + `tests/scripts` **1056 passed, 1 skipped, 0
failed** measured before any re-pin, i.e. every failure this wave produced was in
`tests/scientific`.

 1. `test_wave_p_chemistry_2026_08.py::test_new_families_collect_no_accidental_ph_or_water_activity_correction`
    — SPLIT IN TWO. The pH half is kept and STRENGTHENED (renamed `..._no_accidental_ph_correction_from_their_name`);
    the aw half was asserting the ABSENCE of a correction that has since been deliberately
    ADDED (5 of the 6 Wave P families are net water-releasing), so keeping it would have pinned
    a defect. Replaced by `test_new_families_water_activity_membership_matches_their_measured_stoichiometry`,
    which pins the measured net-water count per family, not the correction value. The original
    docstring is preserved verbatim inside the new one.
 2. `test_wave_p_chemistry_2026_08.py::test_hofmann1998_after_the_refit` — MFT 283.59 → 154.85,
    FFT 297.28 → 267.50, with the full conflict note and the pointer to the misattributed
    conditions.
 3. `test_wave_p_chemistry_2026_08.py::test_the_second_mft_channel_now_contributes_...` — the
    C2+C3 lane's contribution ratio 1.3054 → 1.2060. The lane still contributes, which is the
    property that test owns.
 4. `test_wave_s1_additive_flux_2026_08.py::test_hofmann_sulfur_pair_after_the_additive_propagator`
    — 283.5889/297.2755 → 154.8459/267.4965, with an explicit note that the movement is NOT the
    propagator and that the additive property is asserted separately.
 5. `test_free_aa_quantitative_regression.py` — MFT band (1.00, 1.206, 1.50) → (1.00, 2.209,
    2.75) carrying the Wave S1 pin's RELATIVE span; FFT (1.00, 1.486, 1.95) → (1.00, 1.338,
    **1.76**), i.e. **TIGHTENED** on the half that improved, so improving does not buy slack.
 6. `test_honest_headline_guards.py::test_honest_external_literature_coverage_is_1_of_3_...`
    — RENAMED to `..._is_0_of_3_...`; hits (1,3) → (0,3) and median CI width 0.9463 → 0.7460
    dex, with the per-row cause (Resconi furfural leaving its interval) and the note that a
    narrowing interval losing coverage is the worst of both. Still two-sided and exact.
 7. `test_honest_headline_guards.py::test_pentose_hexose_mft_ordering_is_7_78x_...` — RENAMED
    to `..._is_18_27x_not_the_retired_7_78x_6_15x_3_39x_8_98x_or_15_8x`; ribose 824.7 → 374.0,
    glucose 106.0 → 20.5, ratio 7.78 → 18.27, structural share 3.14 → 4.27, and the
    README/AUDIT doc-token assertion moved 7.78 → 18.27 (both docs updated in the same edit).
    The failure message says in terms that this is not improved discrimination.
 8. NEW `tests/scientific/test_wave_s1b_ph_aw_routing_2026_08.py` (9 tests) — `get_ph_multiplier`
    is reachable from the prediction path and moves `get_rate_constant`; the gate admits only
    the branch point and refuses the ten families the substrings would have claimed; NO family
    receives the pH physics twice (checked over all 29 emitted families); the pyrazine branch
    reaches `Aminoketone_Condensation` at pKa 6.5; dimethylpyrazine RISES with pH; **aw set
    membership re-derived from the enumerated steps and checked against the shipped sets**;
    the furan track is reached and the Wave H case still returns 1.0; the hydrolysis arm; and
    the seven-key dead census.

NOT RE-PINNED, DELIBERATELY: the Hofmann contract (1.45x / 0.09 dex, now failing both by
more), the 3.0x floor in `test_pentose_hexose_sulfur_ordering.py`, the Pratap-Singh tolerance,
the Trikusuma observability factors, `max_fold_error`, the hold-out pre-widening 1/5, the 0/6
predictive headline, the 4/14 and 6/14 aggregates, strict-ready 0/14.

### (l) DATA CURATION — two notes, no values changed

**Trikusuma.** `data/benchmarks/pea_isolate_uht_140C_Trikusuma2019.json` gains a
`content_verification_note` recording Wave S2's retrieval: values CONFIRMED against the OSU
thesis (OhioLINK `osu1531495328317918`, Table 6 p.35 — 781.72 ± 58.59 / 163.16 ± 15.06 /
23.98 ± 0.80 µg/L), the **Processed** column, with Control (331/59/8) and Aged (683/197/22)
named explicitly so the `heated_matrix` state label cannot be silently re-read as either.
Citation year corrected **2019 → 2020** in the benchmark file, `benchmark_intake_registry.json`,
`slr_incorporation_matrix.json` (x2), `computational_priors.json` and `src/benchmark_labels.py`
(both display spellings). **The `benchmark_id` string is UNCHANGED** — it is referenced by the
intake registry, `src/matrix_calibration_registry.py`, six test modules and every generated
artifact; renaming it is an identity change, not a citation fix. The note says so, and says the
year inside the id is wrong.

**Hofmann.** `data/benchmarks/cys_ribose_140C_Hofmann1998.json`'s existing
`content_verification_note` gains a `wave_s2_followup` block: retrieval exhausted (ACS 403,
Unpaywall `oa_locations: []`, OpenAlex/S2 CLOSED, five OA citing papers in full text, Cerny
2015 quoting no number); the conditions appear MISATTRIBUTED (140 °C / 30 min / equimolar is
Mottram & Nobrega 2002's protocol; Cerny 2015 states Hofmann & Schieberle's aqueous model was
145 °C / 20 min at 1:3, and Mottram & Nobrega reported headspace GC-MS, not solution ppb, so
neither is established as the source of 342/200); and the implied **MFT 0.0300 mol % / FFT
0.0175 mol %** on the UNVERIFIED 10 mM basis, with the forward check (1.4 mol % → 15 984 ppb,
46.7x; 0.05 mol % → 571 ppb, 2.85x). **No value, DOI, condition or contract threshold edited.**
[P] — resolving it needs ILL access to the body or a decision to re-express in mol %.

### (m) DOC SYNC

  `docs/validation/directional_accuracy_report.md`: a dated **ADDENDUM** appended, with the
    Wave S2 text left standing verbatim above it — the pre-fix numbers are the evidence. It
    carries the baseline correction, the three confirmations, the disjointness argument, the
    seven-key census, both self-inflicted bugs, the pre/post tables, the per-claim flips, the
    lost row's mechanism and owner decision, the structural reason aw is still 0/3, the four
    degraded benchmarks with per-fix attribution, and a restatement that §7's licence is
    unchanged.
  `README.md`: calibration badge and headline **1/3 → 0/3** with the interval-width table
    (0.95 → 0.75 / 2.21 → 2.96 / 3.56 → 4.18 dex) and the "narrower AND missing" reading; a new
    **Wave S1b** block after the Wave S1 table with the directional pre/post, the four degraded
    rows, the Hofmann conflict, and the structural reason aw is still 0/3; pentose block
    7.78x → **18.27x** with the honest "both numbers fell" explanation.
  `AUDIT.md`: the "what survives is ordering" headline 7.78x → 18.27x with the full split
    history and the denominator-collapse explanation; 13 new rows in the "What Round 3 cost"
    table, every one labelled *(Wave S1b)*, including all four degradations and the 1/3 → 0/3;
    a new section **"Wave S1b — the pH and water-activity physics that was written but never
    connected"**.
  `docs/reference/VALIDATION_CONTRACT.md` §3E: a dated paragraph recording that Wave S1b left
    all eight hold-out points bit-identical AGAIN, with the structural reason and the
    conclusion that two consecutive waves of prediction-path changes have now been invisible to
    the hold-out.
  README forward-mode SAMPLE TABLE re-copied from a live
    `generate_report_visual_examples` run, never by hand, as its own caption requires.
    **2,5-dimethylpyrazine dropped out of the top five and bis(2-methyl-3-furyl) disulfide
    took its place** — at that example's pH 5.5 the pyrazine step is penalised twice (amine
    mostly protonated at pKa 6.5 AND two waters shed), so it falls ~450x. Same swap as the
    snapshots' ranking contract. A one-line note in the Wave S1b block says so.
  CITATION-YEAR PROPAGATION: `src/matrix_calibration_registry.py`'s three
    `source="Trikusuma 2019 ..."` provenance strings and one comment corrected to 2020, then
    `family_deviation_audit`, `validated_envelope_report`, `trace_key_values`,
    `literature_learning_loop`, `literature_backlog`, `matrix_benchmark_deltas` and
    `report_visual_examples` regenerated so no generated artifact still prints the wrong
    year. DELIBERATELY LEFT: `results/validation/citation_verification_ledger.json`, which
    records the earlier finding "CrossRef issue year 2020 vs repo 2019 (online-first) -
    benign" — that is a historical record of the discrepancy being seen and waved through, and
    Wave S1b's correction is more legible standing next to it than replacing it. Also left:
    `results/validation/matrix_observability_refit_pratap_singh.*`, because re-running a refit
    generator to fix a prose year is not a trade this wave will make.

### (n) [P] CARRIED FORWARD

 1. **The projection budget has no moisture dependence**, and HMF shares its producing family
    with furfural, so no family-level aw correction can move HMF against aw in a two-precursor
    system. This is why AW-01/AW-03 are unfixable by routing.
 2. **`Fructofuranosyl_Dehydration`'s pH treatment is an open owner decision** — see (g). Two
    routes to the same product currently get opposite pH treatment.
 3. The aw inhibition applies once per step regardless of how many waters the step sheds; the
    lumped-redox families' Le Chatelier argument is weaker than the pure condensations'.
 4. **The acrylamide destruction term is still missing** (AW-02, TEMP-01, ACR-02). Wave S2
    called it the most consequential single miss; it is still open.
 5. **The external hold-out still cannot see the Maillard network** — two consecutive waves of
    prediction-path change have moved zero of its eight points.
 6. **The Hofmann benchmark's 342 / 200 ppb targets are a repo-internal derivation from an
    abstract-reconstructed range table** (Wave S2b, section below), so the panel's TIGHTEST
    contract is ~1.7x tighter than its own source band and the MFT > FFT ordering it encodes
    is midpoint selection rather than measurement. That is the yardstick against which this
    wave's largest single degradation is scored. Wave S2b's recommendations (mark
    `no_verifiable_source`, retire the 1.45x / 0.09 dex contract, demote out of PRIMARY, flag
    Wave P's refit as anchored on a repo-internal midpoint) are unexecuted OWNER decisions.
 7. `pathway_ranker.py`, `kinetics.py` and `cantera_export.py` still call `get_ph_multiplier`
    directly with the raw substring matching, including the dangerous `"condensation"` key.
    Those paths do not feed `evaluate_benchmark_payload`, so they were out of scope, but they
    are now the only remaining substring consumers of that function.
 8. Everything Wave S1 carried forward that this wave did not touch.

### (o) GATES + SUITE

  GATES, all three re-run on the final tree, after every code, artifact, test and doc edit:
    citation_gate   PASS — 82 files, 964 DOI-bearing fields, 317 unique DOIs (unchanged;
                    this wave introduced NO new citation — the Trikusuma edit corrected a
                    YEAR, not a DOI).
    holdout_guard   PASS — 3/3 invariants, re-run after the regeneration.
    fit_target_gate PASS — both checks. No new fit declared, no fitted constant moved.

  INTERMEDIATE RUNS, kept because they isolate the blast radius:
    - `tests/unit` + `tests/integration` + `tests/scripts` after all three fixes and the
      regeneration but BEFORE any re-pin: **1056 passed, 1 skipped, 0 failed** (525.49 s).
      Every failure this wave produced was in `tests/scientific`.
    - `tests/scientific` at that same point: **10 failed, 208 passed, 2 xfailed** (730.78 s).
      All 10 are the re-pins itemised in (k); none was a code defect.

  FULL SUITE #1 (`tests/unit tests/scientific tests/integration tests/scripts`, documented
  conda path): **1284 passed, 1 skipped, 2 xfailed, 0 FAILED** in 892.75 s, exit code 0.
  SUPERSEDED, and here is why rather than quoting it silently: after it finished, the
  Trikusuma citation year was propagated into `src/matrix_calibration_registry.py` and seven
  more generated artifacts, and the README forward-mode sample table was re-copied from a live
  run.

  FULL SUITE #2: **1284 passed, 1 skipped, 2 xfailed, 0 FAILED** in 1015.92 s, exit code 0,
  reproducing run #1's counts EXACTLY. Also superseded — the Wave S2b reconciliation landed
  after it.

  CERTIFYING FULL SUITE #3, on the final tree with nothing edited afterwards:
  **1284 passed, 1 skipped, 2 xfailed, 0 FAILED** in 1042.96 s, exit code 0, zero FAILED/ERROR
  lines — the same counts a third time, across three trees that differ only in prose, comments
  and one JSON note, which is itself the evidence that none of those edits was load-bearing.
  The 1 skip and both xfails are the declared, strict-marked ones from Wave J2
  (`xfail_strict = true`), so neither can silently start passing. The 1284 = Wave S1's
  certified 1274 + 9 new tests in `tests/scientific/test_wave_s1b_ph_aw_routing_2026_08.py`
  + 1 from splitting the Wave P pH/aw guard into a pH test and a stoichiometry test.
  Arithmetic: Wave S1's certified 1274 + 9 new tests + 1 from the Wave P guard split = 1284.
  All three gates were re-run green on this same final tree after the last edit.
  NOT COMMITTED, NOT STASHED — handed to the orchestrator as instructed.

  RECONCILIATION WITH WAVE S2b, executed after the runs above and re-certified. Wave S2b
  landed in this ledger while Wave S1b was regenerating, and it names Wave S1b's own prose:
  the "read the Hofmann row as a conflict, not a bug" framing rests on a non-measurement.
  VERIFIED INDEPENDENTLY (`data/benchmarks/maillard_validation_benchmarks.md:88` reads
  `| Ribose + Cys, pH 5 aqueous | 140 | 30 min | ~0.02–0.05 | ~0.01–0.03 |`) and CORRECTED in
  five places — README, AUDIT, the report addendum, and the two test comments in
  `test_wave_p_chemistry_2026_08.py` and `test_free_aa_quantitative_regression.py` — plus the
  benchmark file's `wave_s2_followup` note, which now carries a `superseded_by` pointer, Wave
  S2b's repo-internal-origin finding, its off-by-one correction (Cerny 2015 cites Hofmann &
  Schieberle **1995**), and an explicit list of the structural changes Wave S2b recommends
  that NO wave has executed. NO assertion was changed by this reconciliation and no contract
  was relaxed; only prose and comments moved. `tests/scientific` re-run on the corrected tree:
  **228 passed, 2 xfailed, 0 FAILED** (520.89 s); all three gates re-run PASS; then the final
  full suite below.

  ENVIRONMENT NOTE. Free RAM fell to ~75 MB during the long runs and pytest advanced a few
  tests per minute through the `tests/scientific` governance region — Wave S1's and Wave P's
  memory-pressure condition, NOT the ENOSPC condition. It slows the suite; it produced no
  spurious failure in any run this wave.

## Wave S2b — the 342 / 200 attribution, settled (2026-08-27, read-only literature wave)

**No repo file was edited by this wave except this ledger entry.** Full dossier (retrieval log,
re-anchoring analysis, ILL request pack):
`/private/tmp/claude-501/-Users-pabloantoniomorenocasares-Developer-Maillard/dfacf863-78ee-4513-826f-82e7cb2949c6/scratchpad/hofmann_vs_mottram_anchor_dossier.md`

### (a) HEADLINE — the origin of MFT 342 ppb / FFT 200 ppb is INSIDE THIS REPOSITORY

Neither Hofmann & Schieberle 1998 nor Mottram & Nobrega 2002 is the source. The values were
derived from `data/benchmarks/maillard_validation_benchmarks.md` §1.3, an
abstract-reconstructed range table committed in `c7efbbc` — the **same commit** that created
`data/benchmarks/cys_ribose_140C_Hofmann1998.json`. §1.3 tabulates, for "Ribose + Cys, pH 5
aqueous | 140 | 30 min", **MFT `~0.02–0.05` mol % and FFT `~0.01–0.03` mol %**. On the file's
declared 10 mM basis with MW 114.17 (MFT and FFT are both C5H6OS):

* MFT **0.0300 mol %** -> 3.0e-4 x 0.010 M x 114.17 = 3.4251e-4 g/L = **342.5 -> 342 ppb**
* FFT **0.0175 mol %** -> 1.75e-4 x 0.010 M x 114.17 = 1.998e-4 g/L = **199.8 -> 200 ppb**
  (0.017321 = the exact geometric mean of the `~0.01–0.03` band -> 197.7 ppb)

Both recorded values are interior points of the two invented bands. Confidence ~90%.

§1.3 is abstract-derived guesswork, on four independent tells: (1) its only bold, unhedged row
(`1.4` / `0.05`) is verbatim from the abstract; (2) its `Furfural + H2S ... ~0.5` cell is
arithmetic on the abstract's "10 times higher efficiency" (10 x 0.05), but assigned 140 C /
30 min when the abstract's intermediate systems are 180 C / 6 min anhydrous; (3) its
`Glucose + Cys` cell reads `~10x lower than ribose` -- prose in a numeric column; (4) the
140 C / 30 min conditions are transplanted from §1.1 of the same document, which is Mottram &
Nobrega 2002, via §1.1's false premise sentence: *"Fully quantitative SIDA data for MFT/FFT
**from the same system** available in Hofmann & Schieberle (1998)"* (whose cross-reference,
"see §1.4 below", points at Brands & van Boekel -- another generation artifact).

### (b) Three consequences that land on shipped code and shipped prose

1. **The 1.45x / 0.09 dex contract is 1.7x tighter than its own source band.** The md's MFT
   band is 0.02–0.05 mol % = **228–571 ppb, a 2.5x spread**. This disqualifies the contract
   independently of the attribution question.
2. **The README's "read the Hofmann row as a conflict, not a bug" passage rests on a
   non-measurement.** The two bands OVERLAP (MFT 228–571 ppb, FFT 114–342 ppb). The
   MFT 342 > FFT 200 ordering is an artifact of picking 0.030 and 0.0175 inside them; taking
   the FFT band's top and the MFT band's bottom reverses it. The pH-5 conflict may not exist.
3. **Wave P's `thiol_addition_pentodiulose` 28.60 -> 26.35 kcal/mol refit is anchored on a
   repo-internal midpoint**, having been fitted against this benchmark alone. Either revert to
   the un-fitted class value 28.60 or carry an explicit provenance warning until ILL lands.

### (c) Mottram & Nobrega 2002 CANNOT be the absolute anchor -- do not re-anchor to it

Citation verified against CrossRef exactly as recorded (JAFC 2002, 50(14), 4080–4086,
`10.1021/jf0200826`, PMID 12083887). Full text CLOSED on every route (Unpaywall
`oa_locations: []` and `has_repository_copy: false`; OpenAlex/S2 closed; CentAUR does not
resolve locally, was reached via public-DNS `--resolve` to `reading.eprints-hosting.org` and
is behind Anubis bot protection -- and holds no copy anyway; five OA citing papers pulled in
full text quote no value). **Every quantitative statement in its abstract is comparative** --
"much lower quantities", "similar quantities", "relatively unreactive" -- with no unit, no
temperature and no time. It is a headspace-on-Tenax measurement, i.e. a trapped mass, not a
solution concentration; converting one to the other needs a 140 C air/water partition
coefficient, purge efficiency and breakthrough behaviour the paper cannot supply. The repo's
own §1.1 already says *"Concentrations are semi-quantitative (peak areas)."* ~85% confident,
pending the Methods extract.

**What it CAN anchor:** ordinal claims -- ribose-5-P ~= ribose >> IMP, buffered >> unbuffered
for 2,3-enolisation products. Cerny 2015 independently confirms its protocol verbatim:
*"equimolar amounts of L-cysteine and ribose (at 140 C for 30 min)"*. So the JSON's
`conditions` block is a faithful copy of the **wrong paper's** protocol.

### (d) CORRECTION to the `wave_s2_followup` note in the benchmark file

That note says Cerny 2015 attributes 145 C / 20 min at 1:3 to *Hofmann & Schieberle*. **Off by
one paper.** Cerny's sentence cites **Hofmann and Schieberle, 1995** = *Evaluation of the key
odorants in a thermally treated solution of ribose and cysteine by AEDA*, JAFC 43, 2187–2194
(`10.1021/jf00056a042`). His reference list disambiguates: 1998a = the Z. Lebensm. dry-roasting
paper; **1998b** = our `10.1021/jf9705983`. Cerny states nothing about our paper's conditions.
That 10.1021/jf9705983 reused the 145 C / 20 min aqueous protocol is an **inference** (~85%),
supported only by a secondary, search-engine-mediated sentence from BBB 2016
(`10.1080/09168451.2016.1238295`), whose full text tandfonline would not serve.

### (e) RECOMMENDATION -- staged (iii) then (i), with Mottram demoted to ordinal

**Now, without waiting for ILL:** mark 342 / 200 `no_verifiable_source` (repo-internal
derivation, with the §1.3 trail written into the file); **retire the 1.45x / 0.09 dex
contract** and demote the benchmark out of PRIMARY and out of the strict-ready panel; flag the
Wave P refit per (b)(3); correct the README passage per (b)(2).

**After ILL:** rebuild from Hofmann's actual table in **native mol %, not ppb**. mol % is
basis-free; the current ppb target smuggles in the unattested 10 mM basis as a **free
multiplicative parameter sitting underneath the panel's tightest contract**. New
`benchmark_id` reflecting the printed conditions. If the extract shows no aqueous
ribose/cysteine row, fall back to Hofmann & Schieberle 1995; if that also fails, **retire the
absolute sulfur anchor and say so** -- an honestly uncalibrated sulfur branch beats this.

**Mottram & Nobrega -> a separate ORDINAL-only benchmark** (restricted option (iv)): one
absolute anchor (Hofmann, mol %), one ordinal (Mottram, orderings), never two absolute ones.
This is a net gain -- the directional panel is 4/7 on pH and 0/3 on moisture and has no clean
buffer-on/off ordinal test.

### (f) ILL pack for the owner (both primary citations CrossRef-verified, no correction needed)

1. **Hofmann & Schieberle, JAFC 1998, 46(1), 235–241, `10.1021/jf9705983`, PMID 10554225.**
   Photograph: the Experimental section covering the **aqueous** reactions (mmol of cysteine
   and of each carbohydrate, buffer identity/molarity/pH, **total liquid volume**, vessel, T,
   t); **every yields table complete with headers, footnotes and caption** (the header carries
   the unit string; the footnote names the yield's denominator); the ribose/xylose/glucose/
   rhamnose rows and any pH-series rows; the SIDA labelled-standard and response-factor
   statement. **Q1:** is there a quantitative aqueous ribose+cysteine MFT/FFT yield, and at
   exactly what T / t / pH / buffer / cys:ribose ratio (145/20 at 1:3, or 140/30 at 1:1)?
   **Q2:** what is the native unit and denominator, and does the Experimental section give the
   volume and moles needed to convert mol % to a concentration at all?
2. **Mottram & Nobrega, JAFC 2002, 50(14), 4080–4086, `10.1021/jf0200826`, PMID 12083887.**
   Photograph: Methods in full (model-system prep AND the headspace/GC-MS section -- trap
   material, purge flow, collection time, desorption, **internal standard identity and
   amount**); **the sentence defining what the tabulated numbers are** (this is the single most
   important item in the pack); every volatile table complete, with the MFT and FFT rows across
   the full ribose / R5P / IMP x buffered / unbuffered grid. **Q1:** are the numbers ABSOLUTE
   (mass or concentration, defined denominator, response factors applied) or RELATIVE (peak
   area / IS-normalised, no response factor)? **Q2:** exact T, t, amounts, molar ratio, buffer,
   pH, volume.
3. *(only if Item 1 Q1 returns "no aqueous row")* **Hofmann & Schieberle, JAFC 1995, 43(8),
   2187–2194, `10.1021/jf00056a042`** -- Experimental section plus any quantitation table.
   Does it report absolute MFT/FFT concentrations, in what unit?
4. *(optional, if ACS ILL is slow)* **Nobrega & Mottram, Developments in Food Science, 1998,
   483–492, `10.1016/s0167-4501(98)80070-8`** -- the conference precursor to Item 2, same
   experiment, often easier to source as a book scan and likely to carry the same tables and
   method text.

## Wave S2c — the tightest contract was anchored to the repo's own guess (2026-08-27, owner-approved execution)

Executes Wave S2b's staged plan **(iii)**: retire now, rebuild later from the owner's ILL
retrieval. Built on top of Wave S1b's uncommitted tree; every S1b edit is intact and this
wave's changes sit on top of it.

**THE HEADLINE, stated plainly because it is the whole wave:
THE SULFUR BRANCH NOW HAS ZERO ABSOLUTE LITERATURE ANCHORS.**
After Round 2 quarantined three fabricated sulfur benchmarks, the repo said in six places —
README, AUDIT, `src/barrier_constants.py`, two refit records, four test docstrings and
`src/curated_pathways.py` — that one literature anchor survived. It never existed.
`cys_ribose_140C_Hofmann1998`'s MFT 342 ppb and FFT 200 ppb are interior points of two
**invented** mol % bands in `data/benchmarks/maillard_validation_benchmarks.md` §1.3, an
abstract-reconstructed table committed in `c7efbbc` — the *same commit* that created the
benchmark JSON. On the file's declared (and itself unattested) 10 mM basis with MW 114.17:
`0.0300 mol % x 0.010 M x 114.17 = 3.4251e-4 g/L = 342.5 -> 342 ppb`, and the **geometric
mean** of the `~0.01-0.03` FFT band, `0.017321 mol % -> 197.8 -> 200 ppb`. ~90% confidence,
arithmetic exact (Wave S2b).

### (a) THE RETIREMENT RECORDS

`data/benchmarks/cys_ribose_140C_Hofmann1998.json`:

1. **Both `measured_volatiles` values marked `value_status: no_verifiable_source`**, each
   carrying a `content_correction_note` with the §1.3 row quoted verbatim, the exact
   arithmetic above, the four tells that §1.3 is reconstructed, and an explicit
   "do not fit / do not cite / do not report as accuracy against literature".
   **NEITHER VALUE WAS EDITED**, and each says why in a
   `value_left_unchanged_deliberately` field: the shipped number is the evidence that this
   specific number was carried as a literature measurement for four months. The label is what
   this wave corrects, not the value.
2. **The 1.45x / 0.09 dex contract RETIRED.** `validation_contract.scale_thresholds` is
   replaced by a `validation_contract.RETIRED` record holding the retired thresholds, the
   reason (**the contract was ~1.7x tighter than the 2.5x spread of the band its own target was
   interpolated from** — 0.02-0.05 mol % = 228-571 ppb), the overlapping-band finding (MFT
   228-571 vs FFT 114-342, so even the MFT > FFT *ordering* is midpoint selection), and the
   rebuild condition.
   **NOTHING LOOSER WAS INVENTED TO REPLACE IT, AND THE FALLBACK IS STATED RATHER THAN
   HIDDEN:** `_resolve_scale_thresholds` falls back to the global free-precursor defaults
   **1.5x / 0.10 dex**, which are *marginally looser* than what was retired. That is written
   into the file under `WHAT_APPLIES_INSTEAD_read_this_before_calling_it_uncontracted`,
   together with the measurement that it buys nothing: post-revert the row scores
   **4.3797x / 0.4041 dex**, failing the inherited default by ~2.9x on ratio and ~4x on log
   error — a far wider failure than the 2.2086 / 0.2352 it was failing the retired contract by.
3. **`metadata.tier` demoted `PRIMARY` -> `REFERENCE`**, with a `tier_history` block. Tier
   vocabulary was PRIMARY/SECONDARY only; `REFERENCE` is added and documented in
   `src/validation_contract.py` next to `strict_gate_tiers` (which is `("PRIMARY",)`, so the
   demotion is the whole mechanism). `metadata.notes` rewritten: the old text — *"Quantitative
   free amino-acid sulfur benchmark used for strict-ready calibration"* — was false in both
   halves and is quoted in place of itself.
4. `content_verification_note.wave_s2c_retirement` added: what was done, the measured cost,
   what is still true about the file (the DOI is real and CrossRef-verified; the conditions
   block is Mottram & Nobrega 2002's genuine protocol under the wrong citation), a `do_not`
   list, and the rebuild plan. The Wave S2b `owner_decisions_...not_taken` key is rewritten to
   EXECUTED, quoting its own former text.
5. **The file is NOT deleted and NOT moved.** It still executes the whole free-precursor sulfur
   network (thiohemiacetal -> FFT, pentodiulose -> MFT) plus the pH-5 enolisation branch point,
   and it is the only free-precursor cysteine/ribose system in the tree.

`data/benchmarks/maillard_validation_benchmarks.md`:

* **§1.3 heading changed to "⚠️ ABSTRACT-RECONSTRUCTED GUESSWORK — NOT A TRANSCRIPTION"** with
  a dated warning box: the derivation table (mol % chosen -> arithmetic -> shipped ppb -> band
  -> band in ppb), what it caused, the four tells, what the paper's abstract *does* establish,
  and an explicit withdrawal of §1.3's "Benchmark use" instruction (*"Predicted MFT yield …
  must fall in the 0.02-0.05 mol% range"*). **Kept, not deleted — it is the provenance record
  of a fabrication that reached shipped code.**
* **§1.1's false-premise sentence deleted in place and quoted**: *"Fully quantitative SIDA data
  for MFT/FFT **from the same system** available in Hofmann & Schieberle (1998) — see §1.4
  below."* Both halves are wrong: "from the same system" is what licensed transplanting Mottram
  & Nobrega's 140 C / 30 min protocol onto Hofmann's citation (and thence into the JSON, which
  still carries it), and §1.4 is Brands & van Boekel, not Hofmann & Schieberle.

### (b) THE REVERT, WITH THE FULL HISTORY

`src/barrier_constants.py`, `thiol_addition_pentodiulose`: **26.35 -> 28.60 kcal/mol**,
ESTIMATED. The rationale now opens with the ordered history so the middle step cannot
disappear — *(1) 28.60 ESTIMATED (Wave N, the un-fitted `thiol_addition` class value);
(2) 26.35 FITTED (Wave P item 1) against `cys_ribose_140C_Hofmann1998` ONLY; (3) 28.60
REVERTED (Wave S2c)* — followed by the Wave S2b citation and arithmetic, then the Wave P fit's
own text kept **verbatim** under "WHAT THE RETIRED FIT SAID". The Wave K caveat it quoted is
kept and marked SUPERSEDED, with the reason: the question was never whether the mol%->ppb
conversion was documented, it is that **there is no measurement on the far end of it**.
Same treatment as Wave I's Methional revert.

`results/validation/sulfur_barrier_refit_pentodiulose.{json,md}`: **RETRACTED**, mirroring
`hydrolysate_observability_rederivation` exactly — a leading top-level `"RETRACTED"` key in the
JSON (status / retracted_on / retracted_by / reason / consequence / do_not_use_as) and a
`> # ⛔ RETRACTED` blockquote banner above the preserved record in the MD, ending with "the text
below is preserved verbatim as the forensic record". The banner names its own first line as
false.

`results/validation/sulfur_barrier_refit_hofmann.{json,md}`: **ANNOTATED, NOT RETRACTED**, and
the reason is recorded in both files. Its knob `thiol_addition_norfuraneol` = 26.85 sits on a
family **no step emits** since Wave N, so reverting it would move no prediction while making
the provenance record harder to read. Its `contamination_review` verdict *"THIS REFIT IS
UNCONTAMINATED — its result STANDS"* is **withdrawn**: that review checked the fit target was
not one of the quarantined files and never checked whether the target's *values* came from the
paper it cites. Its "standing caveat" about a suspiciously tight contract is marked RESOLVED
and worse than it guessed. `src/barrier_constants.py`'s `thiol_addition_norfuraneol` entry
carries the same annotation: **the Wave H fit was never a literature constraint**, and no
document may describe it as one.

`src/barrier_constants.py` module header: the 2026-08-26 forensics paragraph claiming *"only
thiol_addition retains a literature constraint (Hofmann1998, admissible window [28.10, 28.85]
kcal/mol)"* is kept verbatim under a dated **CORRECTION** — that window is a window around the
repo's own guess.

### (c) EVERY DOWNSTREAM CLAIM CORRECTED (the sweep, itemised)

| File | What it claimed | What it says now |
|---|---|---|
| `README.md` (calibration section) | Wave S2b's retirement recommendation "unexecuted" | new **Wave S2c** block: what was undone (5-row table), the published cost (4-row table), what is *not* fixed |
| `README.md` (Wave P refit passage) | describes the 26.35 refit as shipping | marked **REVERTED by Wave S2c**, read as history |
| `README.md` "The one surviving literature constraint…" | asserted the anchor | struck through, corrected in place: **there is no surviving literature constraint on the sulfur branch** |
| `README.md` projection-refit numbers (1.26x / 0.0497 dex) | live | marked stale on both counts; record deliberately not re-run |
| `README.md` Kinetics section | "Literature-calibrated … anchored to Hofmann, Martins, Nursten" | Hofmann anchor removed; Martins/Nursten stand |
| `README.md` + `AUDIT.md` pentose ≫ hexose | **18.27x** | **8.26x**, with the structural-share reading (23% -> 52%) |
| `AUDIT.md` Round 3 finding #6 | SERIOUS, "reports mol %, not ppb" | escalated to **SERIOUS -> FATAL**: not a units problem, a fabrication |
| `AUDIT.md` "What Round 3 cost" table | ended at Wave S1b | **11 new `*(Wave S2c)*` rows** |
| `AUDIT.md` Open item 1 table | column headed *Measured* | headed *"Measured"*, with a correction box; the structural Wave I finding explicitly still stands |
| `AUDIT.md` | — | new section **"Wave S2c — the tightest contract was anchored to the repo's own guess"** |
| `docs/reference/VALIDATION_CONTRACT.md` strict-ready block | "against a 1.45x contract" | the contract no longer exists; retired-not-widened, tier demoted, still 0/14, fallback stated |
| `docs/reference/VALIDATION_CONTRACT.md` §3E | "two consecutive waves invisible to the hold-out" | **three**, and the sharpest case — a *barrier* revert moved zero hold-out points |
| `docs/validation/directional_accuracy_report.md` | Hofmann row 1.4864 -> 2.2086 as a degradation vs a measurement | addendum: same-day retirement, post-revert 4.3797 / 0.4041; S1b's attribution unchanged |
| `results/validation/content_verification_wave_k.md` item 3 | "tightest contract … on a number nobody has read" | RESOLVED, worse than assumed; ILL ask restated as a *rebuild* ask |
| `scripts/generators/refit_sulfur_barriers_hofmann.py` (docstring + emitted MD header) | "the ONLY surviving literature constraint" | withdrawn in the docstring; the generated MD now emits its own correction blockquote |
| `scripts/generators/refit_thiol_addition_pentodiulose_hofmann.py` | quotes the Wave K conversion caveat | SUPERSEDED block + **DO NOT RE-RUN THIS SCRIPT** against this benchmark |
| `src/curated_pathways.py` | "the only sulfur-branch constant with a surviving literature constraint" | corrected; the **routing decision is explicitly unchanged**, only the claim about the number |
| `tests/unit/test_wave_h_2026_08.py` (module docstring) | "the single surviving literature constraint" | corrected — both of Wave H's calibration acts rested on non-evidence |
| `tests/scientific/test_wave_i_network_chemistry.py` | "the ONLY surviving literature constraint" | corrected; "reference" values relabelled repo-internal |
| `tests/scientific/test_thiamine_fragmentation_benchmarks.py` | "the sole surviving literature constraint" | corrected, plus the compounding note (Cerny's own 2.47 ppb is *also* unverified) |
| `tests/scientific/test_free_aa_quantitative_regression.py` | "the only surviving literature constraint on that branch" | corrected in the module header |
| `src/validation_contract.py` | tier vocabulary undocumented | PRIMARY / SECONDARY / **REFERENCE** documented next to `strict_gate_tiers`, with what the demotion does and does not do |

**ONE DOWNSTREAM CLAIM LEFT UNCORRECTED, AND IT IS UNCORRECTED ON PURPOSE.**
`docs/validation/isotope_topology_evidence.md` (~line 102) still reads *"which Wave I then
recorded as **the only surviving literature anchor on the entire sulfur branch**"*. That file
was on this wave's explicit do-not-touch list, so it was not edited. **It is now the last place
in the repository that still asserts the retired claim** and it needs a one-paragraph correction
from whoever owns that file. Flagged for the orchestrator.

**`fit_target_gate` registration — checked against the gate's actual semantics, and NOT
weakened.** The gate skips records carrying a `RETRACTED` key (`RETRACTION_KEYS`). Before this
wave `cys_ribose_140C_Hofmann1998` was `fit_target_of` **three** records; retracting the
pentodiulose one leaves **two live**: `projection_constant_refit.json` and
`sulfur_barrier_refit_hofmann.json`. So the benchmark **remains a declared fit target**, remains
`fitted_row: true` in `prediction_uncertainty.json`, and remains **out of the honest
literature-coverage numerator *and* denominator. That is the correct outcome and it is why
`sulfur_barrier_refit_hofmann` was annotated rather than retracted:** retracting it too would
have removed the last live registration and let a benchmark with no verifiable source back into
the literature-coverage denominator *as if it were literature*. `honest_literature_coverage`
is **0/3 before and after**, `excluded_fitted_rows` **2 before and after**. Verified by running
the gate: PASS on both checks.

### (d) ARTIFACTS — REGENERATED, AND DELIBERATELY NOT

The tier change alone moves almost nothing (strict-ready was already 0/14). **The revert is
what moved things: it is a change to a shipped barrier, so every prediction-consuming artifact
was stale.** That is the criterion used, stated so the choices can be checked.

REGENERATED (11):
`benchmark_summary.{md,json}` · `validation_overview.{md,json,png}` (+ `docs/assets/`) ·
`family_validation_overview.{md,json}` + `family_benchmark_accuracy.png` / `family_parity.png` ·
`prediction_uncertainty.{md,json}` (n=200, seed 0) · `external_validation_report.{md,json}` ·
`experiment_value_ranking.{md,json}` · `loo_leverage.{md,json}` ·
`gap_heatmap.png` + `experiment_brief_cards.html` · `family_deviation_audit.{md,json}` ·
`matrix_sigma_residual_derivation.{md,json}` · `report_visual_examples/` + the three
`docs/assets/report_*.png`, and the README forward-mode sample table **re-copied from that live
run, not written by hand**, as its own caption requires.

**`prediction_uncertainty` WAS regenerated, against the brief's default.** The brief said not
to unless panel membership changed — it did not (11 benchmarks, 35 matched rows, unchanged) —
but the barrier revert changes the *predictions inside it*, and shipping an MC artifact whose
point predictions no longer match the model is exactly the drift this campaign removes. It cost
4m15s. The headline numbers are **unchanged**: coverage 28/35 (80%), `honest_literature_coverage`
0/3, `excluded_fitted_rows` 2/2. Only interval widths moved (fitted-row median 2.9573 -> 2.9535
dex, internal-synthetic 4.1780 -> 4.1084).

**Bit-identical after regeneration, and reported as such:** `external_validation_report`
(all eight hold-out points; **third consecutive wave**) and `matrix_sigma_residual_derivation`.

**The four internal reproducibility snapshots recovered exactly after the documented refresh**
(`refresh_internal_reproducibility_snapshots.py`, `GENERATOR_TAG` **v10 -> v11** with a leading
Wave S2c block). Before the refresh all four *failed* at 2.403x, which is the drift detector
working. Movement: **MFT x0.416079**, **disulfide x0.420049** (MFT-derived, so it follows), and
**every other propagated row x1.010812** — FFT, furfural and DMP alike. One shared factor across
three chemically unrelated compounds is the signature of the fixed volatile budget
re-normalising, not of three coincidences. Hexanal / Nonanal **x1.000000** (injected by the
lipid lane, never propagated). **The ranking contract did NOT change this time** — furfural >
FFT > MFT > disulfide > DMP holds before and after.

DELIBERATELY NOT REGENERATED, with reasons:
* `projection_constant_refit.{json,md}` and `matrix_observability_refit_pratap_singh.{json,md}`
  — **refit generators**. Re-running one to refresh a number nobody may act on is a
  recalibration event dressed as bookkeeping; Wave S1b made the same call. The README numbers
  they feed are marked stale in place, and the decision they record ("do not move the global
  scale") is unchanged and now better supported.
* `sulfur_barrier_refit_hofmann.{json,md}` — annotated, not re-run; its family is dead and
  re-running it would fit against the same non-measurement.
* `citation_verification_ledger.*` — a historical record of what was seen and waved through.
* `content_verification_wave_k.md` — a report, annotated in place.

### (e) TESTS RE-PINNED — 13 changes, every one with a dated causal comment, none relaxed

1. `tests/unit/test_wave_h_2026_08.py::test_corrected_mft_route_barriers_are_estimates_not_fits`
   — 26.35 -> **28.60**; the test returns to its original meaning (ESTIMATE, not fit). Now also
   asserts the rationale carries the **whole** history (`ESTIMATED` / `FITTED … Wave P item 1` /
   `REVERTED … Wave S2c`), `ZERO ABSOLUTE LITERATURE ANCHORS`, and `RETRACTED`, with an explicit
   `!= 26.35` guard whose message says why refitting there is circular.
2. `..._wave_h_2026_08.py` module docstring — the "single surviving literature constraint"
   correction.
3. `tests/scientific/test_wave_p_chemistry_2026_08.py::test_pentodiulose_barrier_is_the_wave_p_fit_and_carries_the_conversion_caveat`
   — 26.35 -> 28.60; keeps every Wave K caveat assertion (they must survive verbatim) and adds
   `SUPERSEDED` + `RETRACTED`.
4. `…::test_hofmann1998_after_the_refit` — MFT 154.85 -> **78.0867**, FFT 267.50 -> **293.6735**.
5. `…::test_the_second_mft_channel_now_contributes…` — 128.40 -> **46.82**, 154.85 -> **78.09**,
   and the C2+C3 lane's contribution ratio **1.2060 -> 1.6678**: the lane got *more* important,
   because the reverted barrier sits on the pentodiulose lane only.
6. `tests/scientific/test_wave_s1_additive_flux_2026_08.py::test_hofmann_sulfur_pair_after_the_additive_propagator`
   — 154.8459/267.4965 -> **78.0867/293.6735**, with the third distinct cause spelled out.
7. `…::test_both_hofmann_mft_routes_share_the_amadori_rate_limiting_step` — **RENAMED** to
   `test_the_hofmann_lane_rate_limiting_steps_are_measured_not_assumed`, because **they no
   longer share it**: see (f).
8. `tests/scientific/test_free_aa_quantitative_regression.py` — MFT band
   (1.00, 2.209, 2.75) -> **(1.00, 4.380, 5.45)**, FFT (1.00, 1.338, 1.76) ->
   **(1.00, 1.468, 1.93)**. Both carry the Wave S1b pins' **relative** spans, so both are
   re-centrings on worse numbers, not loosenings. Header note added: these bands are drift
   guards against a number this repository invented and cannot speak to accuracy.
9. `tests/scientific/test_honest_headline_guards.py::test_pentose_hexose_mft_ordering_…`
   — **RENAMED** `…_is_8_26x_not_the_retired_18_27x_…`; ribose 374.0 -> **169.1**, glucose
   **unchanged at 20.5** (which is what identifies the cause), ratio 18.27 -> **8.2607**;
   README/AUDIT doc-token assertion moved 18.27 -> **8.26** (both docs updated in the same
   edit). The failure message now also says: if it ever *rises*, check first whether something
   was refitted against `cys_ribose_140C_Hofmann1998`.
10. `tests/scientific/test_wave_i_network_chemistry.py::test_hofmann1998_sulfur_predictions_are_recorded_not_tuned`
    — MFT floor 100.0 -> **40.0**, FFT ceiling 500.0 -> **600.0**. **This is a deliberate
    loosening of a drift guard and it is labelled as one**: these bounds never asserted
    accuracy (the accuracy pin is (8), re-pinned *worse*), the new floor keeps its original
    character at roughly half the current value, and the docstring says all of it.
11. `tests/scientific/test_thiamine_fragmentation_benchmarks.py` — Cerny 23.406x -> **25.741x**
    recorded; **the [12, 46] band is deliberately NOT widened** (25.741 sits comfortably
    inside, so the guard keeps exactly its sensitivity). Plus the "sole surviving literature
    constraint" correction and the compounding note.
12. `tests/scientific/test_validation_contract.py` — two tests **retargeted** off the demoted
    file: the strict-gate fixture to `acrylamide_spi_extrusion_130C_ACSRef3` (PRIMARY,
    free_precursor), **and the loss turned into a new guard** — the test now asserts
    `tier == "REFERENCE"` and `strict_ready is False` on Hofmann, so re-promoting it fails
    here. The threshold-override test moved to
    `furosine_extrusion_crossover_140C_RamirezJimenez2000` **and its direction flipped**
    (loosen-then-pass at 1.70, plus a reject case at 2.50), because of a finding worth
    recording: **after the retirement, no free-precursor benchmark in the tree carries a
    contract tighter than the global default any more** (acrylamide 1.5/0.20, cml_cel 1.8/0.25,
    furosine 2.0/0.30, Bolton 3.0/0.48). It also asserts the retired contract **stays** retired
    (no `scale_thresholds` key; the retired values preserved).
13. `tests/scientific/test_benchmark_summary.py` — two fixtures retargeted off the demoted file
    for the same mechanical reason; they were using it only as a generic PRIMARY host and their
    subject is unchanged.

### (f) ONE STRUCTURAL FINDING THIS WAVE PRODUCED, FOUND BY A TEST RATHER THAN LOOKED FOR

**The MFT lane no longer shares its rate-limiting step with the FFT lane.** At 28.60 kcal/mol,
plus the ~0.935 kcal/mol water-activity term at this benchmark's aw 0.98,
`Thiol_Addition_Pentodiulose` has an **effective barrier of 29.5354** — above the shared
`Amadori_Rearrangement` trunk at **29.0603**. Measured on the shipped tree:

```
MFT: Schiff 21.6070 | Amadori 29.0603 | Enolisation_2_3_Amadori 28.0311 |
     Deoxyosone_Reduction 28.9354 | Thiol_Addition_Pentodiulose 29.5354  <- slowest
FFT: Schiff 21.6070 | Amadori 29.0603  <- slowest | Enolisation_Intermediate 21.0000 |
     Enolisation_1_2 27.6964 | Thiohemiacetal_Formation 23.3000 | Thiol_Dehydration 27.7354
```

**It does NOT reopen the propagator rule** — the shipped channel id is the route's full ordered
step-set, so routes are distinct whether or not they share a bottleneck, and the additive sum is
unaffected. **It DOES retire Wave S1's argument** that "Wave P's 242.38 + 71.02 arithmetic never
applied because both lanes bottleneck on the same trunk". That is now half true, and the
surviving half is the decisive one: **the volatile budget is fixed**, so single-lane figures
were never addable wherever the bottleneck sits. The budget invariant is asserted independently
(`test_the_volatile_budget_still_caps_the_sum`, 1 part in 1e12).

### (g) MEASURED MOVEMENT — the cost, published rather than absorbed

| | pre (Wave S1b) | **post (Wave S2c)** |
|---|---|---|
| Hofmann MFT (vs the repo-internal 342) | 154.85 ppb, 2.2086x under | **78.0867 ppb, 4.3797x under** |
| Hofmann FFT (vs the repo-internal 200) | 267.50 ppb, 1.3375x over | **293.6735 ppb, 1.4684x over** |
| Hofmann max_ratio / MALE | 1.45x / 0.09 dex contract, failing at 2.2086 / 0.2352 | **contract retired**; failing the inherited 1.5 / 0.10 at **4.3797 / 0.4041** |
| Cerny 2008 MFT (reference-only) | 23.406x | **25.741x** |
| Pentose >> hexose ordering | 18.27x (4.27x structural, 23%) | **8.2607x** (4.27x structural, **52%**) |
| Worst reference-only point | 23.406x | **25.741x** |
| Median authoritative compound max-ratio | 112.227 | **113.394** |
| SLR-01 family mean abs log10 | 0.283 | **0.456** |
| SLR-09 family mean abs log10 | 0.246 | **0.249** |
| Strict-ready | 0/14 | **0/14** |
| `status_counts` | scale-gap 8 / pass-no-ranking 2 / pass 4 | **unchanged** |
| MC coverage | 28/35 | **28/35** |
| `honest_literature_coverage` | 0/3 | **0/3** |
| External hold-out (8 points) | median 93.68x, coverage 3/8 | **bit-identical** |
| Worst experimental point (Bolton MFT) | 6730.854x | **unchanged** (glucose lane, no pentodiulose step) |

**The one thing that got better is the one that matters.** Only **4.27x** of the pentose >>
hexose ordering was ever structural and that number did not move, so the share of the claim
carried by mechanism rather than by a barrier gap went **23% -> 52%**, and what is left riding
on a gap is **1.05 kcal/mol between an estimated class value and an unconstrained legacy fit**,
where it used to be 3.30 kcal/mol between a **fitted** barrier and that same legacy fit. No part
of the surviving ordering claim now traces to a number this repository invented. The one Hofmann
claim that survives is *pentose >> hexose* itself, because it is stated in the paper's
**abstract**, which is retrievable — unlike the yields table, which is not.

### (h) TESTING — targeted subsets only, per the owner directive. NO full-suite run.

| Tier | Result |
|---|---|
| `tests/scientific` (whole tier — it is the blast radius) | **226 passed, 2 xfailed, 0 failed** (514.96 s) after the re-pins |
| `tests/unit` | **926 passed, 1 skipped** |
| `tests/scripts` + `tests/integration` | **126 passed** |

Entry-state failures, all itemised in (e) and none a code defect: 7 in `tests/scientific`
(Wave P x3, Wave S1 x2, free-AA regression, Wave I bounds) + 1 honest-headline guard + 2
`test_validation_contract` + 2 `test_benchmark_summary` = **12**, plus the Wave H unit test.

**ENVIRONMENT NOTE, so a re-run is not mis-read.** `tests/unit/test_cli_scripts.py` (7 tests)
fails with `FileNotFoundError: 'python'` unless the conda env's `bin` is on `PATH` — it shells
out to a bare `python`. Confirmed environmental, **not caused by this wave**: with
`PATH=<env>/bin:$PATH` all 7 pass. Interpreter used throughout:
`/opt/homebrew/Caskroom/miniforge/base/envs/maillard/bin/python`.

### (i) GATES — all three re-run on the final tree, after every code, artifact, test and doc edit

```
citation_gate   PASS — 83 files, 972 DOI-bearing fields, 318 unique DOIs immediately after
                this wave's last edit (82/964/317 before it: the +1 file and +1 field are the
                real, CrossRef-verified 10.1021/jf9705983 now cited in
                src/validation_contract.py's tier comment). Re-run later on the shared tree:
                84/973/318 — the extra file/field is Wave T3's, not this wave's.
holdout_guard   PASS — 3/3 invariants
fit_target_gate PASS — both checks. One fit record retracted, one constant reverted; the
                benchmark stays a declared fit target of two live records, so `fitted_row`
                and the coverage accounting are unchanged. The gate is not weakened.
```

**CONCURRENCY, reported because it affects how the numbers above should be read.** A
**Wave T3** agent was editing the same working tree during this wave (`data/lit/`,
`src/headspace.py`, and a file-level warning header prepended to
`data/benchmarks/maillard_validation_benchmarks.md` — which cites this wave's finding and does
not conflict with the §1.1 / §1.3 edits made here; both survive). Consequences for the runs
above, checked rather than assumed:

* A later `tests/scientific` sweep showed **1 failed, 227 passed, 2 xfailed** —
  `test_slr_reference_payloads.py::test_slr_incorporation_matrix_covers_new_track1_sources`,
  failing on `jafc_3c05991_hexanal_binding`'s `incorporation_status`
  `encoded_modeled_shown -> encoded_modeled_shown_pending_reanchor`. That record is Wave T3's
  (`git diff` attributes it to a Wave T3 withdrawal note); **this wave never touched
  `data/lit/`.**
* A later run of the Wave S2c touched-test set showed the two
  `test_honest_headline_guards.py` **no_verifiable_source census** tests failing at
  119 records / 101 in `data/lit` against the pinned 102 / 84. Also Wave T3's: the census
  (`_no_verifiable_source_census`) counts objects carrying **`source_status`** under
  `data/lit/**` and `data/qm/**` ONLY. This wave's new markers are **`value_status`** fields
  inside `data/benchmarks/cys_ribose_140C_Hofmann1998.json`, which the census does not scan —
  and that file contains **zero** `source_status` keys. Verified by grep, not by inference.
  Both tests **passed on this wave's own tree** before Wave T3's `data/lit` edits landed.
* **This wave's own touched-test set is green:** the 11 files listed in (e) run
  **81 passed, 0 failed**.

NOT COMMITTED, NOT STASHED — handed to the orchestrator, on top of Wave S1b's intact
uncommitted work.

### (j) [P] CARRIED FORWARD

1. **The benchmark is still SCORED.** It keeps `measured_volatiles`, so its rows are still
   enumerated in `benchmark_summary`, `validation_overview` and `prediction_uncertainty`, where
   they read as large misses **against a fabricated yardstick**. The tier change removes
   strict-gate eligibility and nothing else, and the file says so at length in
   `metadata.tier_history`. Taking it out of the scored population needs an owner decision:
   move the block to `reference_volatiles` (the treatment
   `thiamine_cys_xylose_145C_Cerny2008` already carries for its own unverified value), or
   quarantine the file. **Both change panel membership**, which is why this wave took neither.
2. **The rebuild is the real fix, and it waits on ILL.** Wave S2b §(f) has the pack. Target:
   the paper's own yields table in **native mol %** (basis-free; the ppb target smuggles in the
   unattested 10 mM basis as a free multiplicative parameter under whatever contract sits on
   top), a new `benchmark_id` reflecting the printed conditions, and a contract from the
   paper's own reported precision. If the paper has no aqueous ribose/cysteine row, fall back
   to Hofmann & Schieberle **1995** (`10.1021/jf00056a042`); if that fails too, **retire the
   absolute sulfur anchor entirely and say so**.
3. **Mottram & Nobrega 2002 as an ORDINAL-only benchmark** (ribose-5-P ~= ribose >> IMP;
   buffered >> unbuffered) feeding the directional panel — a net gain, since that panel is 4/7
   on pH and 0/3 on moisture and has no clean buffer-on/off ordinal test. Never as a second
   absolute anchor: it is a headspace-on-Tenax measurement, and the repo's own §1.1 calls it
   semi-quantitative peak areas.
4. **`docs/validation/isotope_topology_evidence.md` still asserts the retired claim** — see the
   note at the end of (c). Do-not-touch list; needs its owner.
5. **The conditions block is still misattributed.** 140 C / 30 min / equimolar / pH 5 is
   Mottram & Nobrega 2002's protocol (Cerny 2015 confirms it verbatim) sitting under Hofmann &
   Schieberle's DOI. Not corrected here, because correcting it means choosing between two
   papers' protocols and that is what the ILL retrieval is for.
6. **Nothing external can see any sulfur constant.** Three consecutive waves — a propagator
   change, a routing repair, and now a *barrier revert* — have each moved zero of the eight
   hold-out points. Combined with (a), there is currently **no evidence of any kind**, in-panel
   or external, bearing on the sulfur-branch barriers.
7. Everything Wave S1b and Wave S2b carried forward that this wave did not touch.

---

## Wave T3 — disarming the laundering engine, and six labelling repairs (2026-08-27)

Specification: Wave T1's forensic sweep (findings T1-00 through T1-09). **No numeric value in
this repository was changed, refitted, substituted or invented by this wave.** Every edit is a
label, a warning, a corrected derivation record, a propagated-and-verified DOI, or a test
re-pin caused by one of those. Two things that were *not* numbers did change: one
`incorporation_status` string and one `provenance_tier`, both demotions, both because a
retrieval settled a question against the repo.

### (a) T1-00 — `scripts/trace_key_values.py` DISARMED, and its report retitled

The single highest-leverage item. The script marked a registry entry **"Fully Verified"** when
every one of its numbers appeared as a *substring* within ±30 lines of a citation *surname*
inside the LLM-generated markdown under `data/Gemini_Deep_Research/`. It never opens a paper,
never contacts CrossRef, and its surname match is unanchored (`Bi` matches `binding`) while its
number match is a bare substring (`4.5` matches `14.52`). Its output published
**"Fully Verified (All values matched): 153 (57.5%)"** under the title *"Numeric Value
Traceability and **Verification** Report"*. Under this repo's own rule
(`data/Gemini_Deep_Research/README.md`: *"the deep-research report says so" is not
provenance*), those 153 rows are the **laundering census**, published as its opposite. It is
the one document in the tree that would have talked a reviewer out of this entire finding.

**OLD → NEW vocabulary, one for one** (the computation is byte-for-byte unchanged; only the
names and labels moved, and the module docstring says so explicitly):

| OLD | NEW |
|---|---|
| `"Fully Verified (All values matched)"` | `"DIGEST-ECHO (NOT VERIFICATION) — every number echoes the LLM digest corpus"` |
| `"Partially Verified (Some values matched)"` | `"PARTIAL DIGEST-ECHO (NOT VERIFICATION) — some numbers echo the corpus"` |
| `"Unverified (No values matched)"` | `"NO DIGEST ECHO — origin unaccounted for even by the digests"` |
| section `"## Fully Verified Entries"` | `"## DIGEST-ECHO entries — every number found only inside LLM-generated text"` |
| `"...are 100% matched and verified in the Deep Research markdown files"` | `"...No paper was opened. This is the list of entries whose only demonstrated upstream is a machine-generated digest — i.e. the remediation worklist, not the safe list."` |
| section `"## Unverified & Partially Verified Entries (Action Required)"` | `"## Entries with NO / PARTIAL digest echo"`, with the note that this says nothing about correctness |
| title `"# Numeric Value Traceability and Verification Report"` | `"# LLM-Digest Echo Census — **NOT a verification result**"` |
| internal `fully_verified` / `partially_verified` / `unverified` | `full_digest_echo` / `partial_digest_echo` / `no_digest_echo` |
| stdout `"Verified: 153/266"` | `"NOTE: this script verifies nothing..."` + `"DIGEST-ECHO (not verified): 153/266"` |

Also added: a ~60-line module docstring stating exactly what the script measures (proximity of
a number to a surname inside LLM-generated text) and what it does not (contact with any paper,
any index, any resolver), plus the sentence a reader most needs — **a high echo count is a bad
sign, not a good one** — and the note that the NO-DIGEST-ECHO class is not thereby better, it
merely has no identified upstream at all. The regenerated report carries the same warning as a
blockquoted header before any number, quoting the Gemini README rule and recording the old
wording so the change is auditable from the artifact itself.

**Counts are unchanged** (266 / 153 / 59 / 54) and deliberately so: re-labelling is not
re-measuring. `results/validation/key_value_trace_report.md` is regenerated on disk; note it is
**gitignored** (`.gitignore:125`, `results/validation/*`), so the tracked artifact of this fix
is the script.

One incidental bug fixed while there: `ROOT = Path(".").resolve()` was the process CWD, so the
script silently produced an empty census unless run from the repo root. Anchored to `__file__`;
verified by running it from `/tmp` and getting identical counts.

**Grep sweep for the 153 / 57.5% figure being cited as evidence: CLEAN.** `README.md`,
`AUDIT.md`, `docs/reference/VALIDATION_CONTRACT.md`, the ledgers and the tests contain no
reference to it. The only mention anywhere was `tasks/audit_remediation.md:4505`, naming
`trace_key_values` in a list of artifacts Wave S1b regenerated — not a claim, left alone.

**Citation gate — check 6 added, and the baseline stayed EMPTY.** `scripts/ci/citation_gate.py`
gains **digest-as-provenance**: a record may not name a machine-generated or
abstract-reconstructed in-repo document (`data/Gemini_Deep_Research/**`,
`docs/research/archives/`, `raw/NN_*.md`, "Literature Report N",
`maillard_validation_benchmarks.md`, `Maillard_meat.md` / `Maillard_Plant_based.md`, the Elicit
dumps, the nonexistent `pathways.md`) in its `citation` / `source` / `source_citation` while
claiming a `provenance_tier` that asserts a primary source was read (`direct_measurement`,
`literature_derived_transfer`, `literature_bounded_provisional`). Tiers that already admit no
paper stands behind the number — `repo_literature_synthesis`, `mechanistic_surrogate`,
`direct_model_assumption`, `unsourced_withdrawn` — are **exempt by design**: naming the digest
there is disclosure, not laundering (this is why `computational_priors.json`
`matrix_corrections[4]`, which cites `docs/Maillard_Plant_based.md` under
`repo_literature_synthesis`, correctly does not fire). `source_status: no_verifiable_source`
also satisfies the check, exactly as in check 5. Separately, `digest_echo` / `llm_digest` /
`deep_research_says` are refused outright as `provenance_tier` or `source_status` values.
This fits the gate's existing semantics (structural, offline, disclosure-enforcing) and
weakens nothing.

It fired on **exactly one** record in the whole tree, and that record was **fixed at source
rather than waived** — see (g) — so `WAIVERS` remains `()` and
`tests/unit/test_audit_remediation_carried_2026_08.py::test_citation_gate_baseline_is_empty`
passes unmodified.

### (b) T1-01 — the runtime MOCK is now labelled and surfaced

`data/lit/protein_source_registry.json` has always described itself as *"Mocked values for 14
protein sources based on Report 06 requirements. Pending empirical substitution."* — where
"Report 06" is `data/Gemini_Deep_Research/06_alternative_proteins.md`. That sentence sat in a
JSON field nothing read, while the numbers underneath it went live:
`src/literature_runtime.py:50` and `src/matrix_correction.py:49` load the file, and
`hydrolysate_observability_bias`, `off_note_penalty` and `lox_activity_flag` enter
`matrix_uncertainty_factor` directly, `meaty_potential_multiplier` drives the meaty-potential
score, and `get_protein_source_profile` reaches `src/recommend.py`.

* **File:** file-level `source_status: no_verifiable_source`, `value_basis:
  mocked_placeholder`, a `provenance` line naming the digest, a `provenance_note`, and a
  ~200-word `_WARNING` field stating that the plant-source DIFFERENTIATION this file encodes
  (pea more off-note-prone than soy, mycoprotein meatier than oat, …) **is not evidence**.
  Every one of the 14 profiles gets the same two labels. **All 14 × 5 numeric values verified
  bit-identical before and after** by an automated diff in the writing script.
* **Runtime, following the family-12 unit-warning precedent** (`_resolve_concentration_unit`):
  both consumers now emit a `RuntimeWarning` at load naming the defect and what depends on it,
  and rescaling nothing. `src.matrix_correction._warn_if_registry_unsourced` is written as a
  reusable predicate returning whether the payload admits it has no source.
* **Surfaced in the payload, not only on stderr:** the family-06 lane payload gains
  `protein_source_provenance` (registry path, status, declared upstream, the list of outputs it
  contaminates, the warning) and `protein_source_profile_unsourced`, so a consumer of
  `matrix_uncertainty_factor` cannot read source differentiation as evidence.
* **NOT DONE, on purpose:** no replacement value invented; `matrix_uncertainty_factor` not
  widened for mock-sourced profiles. Widening it is a science decision — **[P], owner**.

### (c) T1-02 — §3.2's `log_slope = 0.235`: value KEPT, recorded derivation CORRECTED

The record said the surrogate was *"anchored to the Karolkowski 2021 **qualitative** pH trend"*
and its `numeric_reference.units` said `qualitative_release_order` — i.e. the file asserted that
no number had been taken from anywhere. That is not what the constant is. `log_slope = 0.235`
is **exactly ln(1.60)/2 = 0.2350018**, the value that makes
`exp(s·(6.0−4.5)) / exp(s·(6.0−6.5)) = 1.42266/0.88913 = 1.60000` — a +60% ratio. `max_factor =
1.6` is the same number again. The 1.60 is the arithmetic midpoint of the *"~55–65% higher"*
band in `data/benchmarks/maillard_validation_benchmarks.md` §3.2, a fabricated table whose
citation is a self-declaring placeholder DOI and whose own hexanal row (340/205 ppb) implies
+65.9%, not the +59% it is labelled with.

**DOI CrossRef-verified this wave before writing it** (`api.crossref.org/works`): the real
paper behind the one real number in §3.2 is **Fischer, E.; Cachon, R.; Cayot, N. "Effects of
extraction pH on the volatile compounds from pea protein isolate: Semi-Quantification method
using HS-SPME-GC-MS." *Food Research International* 150, 110760 (2021),
`10.1016/j.foodres.2021.110760`** — every field (title, three authors, journal, volume, article
number, year) matches what the repo's own 2026-08-26 `doi_repair` block recorded, and the
abstract's closing sentence is verbatim *"hexanal release was found 59% higher with extraction
using pH 4.5 than with pH 6.5"* (confirmed independently via Europe PMC, PMID 34865778).

Changes, values untouched:
* `runtime_surrogate.surrogate_basis` on **both** blocks (`aldehydes/pea_protein[1]`,
  `furans/pea_protein[0]`) rewritten to state the actual construction, its fabricated origin,
  and the independent vindication — plus `log_slope_derivation: "ln(1.60)/2 = 0.2350021"` and
  `source_status: no_verifiable_source` on the surrogate blocks.
* **The vindication is stated honestly rather than flatteringly.** 1.600 vs 1.590 is 0.63%, but
  Fischer varied the pH at which the isolate was **extracted** (a manufacturing variable
  determining which volatiles the powder carries) whereas `src/headspace.py` varies pH as a
  **release** parameter at measurement time. Different physical quantities; the agreement is a
  coincidence of magnitude, not a measurement of this knob. Both the data record and the code
  comment say so.
* `numeric_reference` on the **hexanal** record corrected from `qualitative_release_order` to
  the real quantity — units, value 59.0, the verbatim quote, where it appears, when it was
  verified, the extraction-vs-release caveat, and a `supersedes` field recording the old
  content. On the **2-pentylfuran** record `qualitative` is *correct* and was kept, with a
  `compound_caveat` recording that Fischer's 59% is hexanal-only ("for example") and that the
  shared `log_slope` was **not** derived from any 2-pentylfuran measurement.
* `src/headspace.py:239` gains a 19-line comment stating the same thing at the point of use.
  The placeholder DOI is quoted there **with spaces inserted** so it is not re-asserted as an
  anchor — the repo's own idiom, and necessary: writing it verbatim made
  `citation_gate` check 2 fire on the comment, which is the gate working correctly.
* **Deliberately NOT re-derived as ln(1.59)/2 = 0.23187.** It would move predictions by <1%,
  it would replace a brief-midpoint with a number measuring a different quantity, and it would
  break the behavioural pin in `tests/scientific/test_matrix_headspace_ph_validation.py` for no
  epistemic gain.

### (d) T1-03 — the compound conflict: RETRIEVED AND SETTLED, and BOTH sides were wrong

`jafc_3c05991_ppi_hexanal_binding` shipped **52.76 ± 4.65 %bound attributed to HEXANAL in PPI
at `provenance_tier: direct_measurement`**. The internal brief's §2.4 attributed ~52.8 ± 4.6 to
**trans-2-heptenal**. This wave retrieved the full text (open access CC-BY, Europe PMC
`PMC10739987` fullTextXML) and **neither is right. It is OCTANAL.** Verbatim, Results:

> "the displacement of the carbonyl group from the inner part of the molecule toward the edge
> leads to a significant binding increase from 14.73 to **52.76% ± 4.65** for the studies with
> plant protein isolates (PPI, SPI, and LPI)"

So 14.73% = 2-octanone, 52.76% = octanal, and it is a **panel-wide range endpoint pooled over
PPI/SPI/LPI**, not a PPI-specific value. The abstract repeats it in the same pooled form. The
brief is wrong too: trans-2-heptenal is **54.60%** (±3.95 / ±2.98 in two comparisons), and its
companion claim that heptanal is ~30–35% is contradicted — heptanal is **13.73%**.

**What the paper says about hexanal: no number, anywhere in the text.** Hexanal appears only in
Figure 4 (storage time × 5/70/90 °C); Tables 1–2 are surface hydrophobicity and sulfhydryl
content, not binding. Its only hexanal statements are qualitative and point the **opposite** way
from how this record was used: *"SPI showed the highest binding capacity for hexanal followed by
LPI … In contrast, PPI and WPI showed the lowest binding affinity for hexanal."*

Correct citation, CrossRef- and Europe PMC-verified: **Barallat-Pérez, C.; Janssen, H.-G.;
Martins, S.; Fogliano, V.; Oliviero, T. (2023). "Unraveling the Role of Flavor Structure and
Physicochemical Properties in the Binding Phenomenon with Commercial Food Protein Isolates."
*J. Agric. Food Chem.* 71(50), 20274–20284.** PMID 38059380.

Actions:
* `provenance_tier` **demoted `direct_measurement` → `unsourced_withdrawn`**, with
  `provenance_tier_history`. It cannot be a direct measurement of a quantity the paper does not
  report.
* `numeric_reference.value` → `null`, `value_status: withdrawn_wrong_compound`, with
  `withdrawn_value` / `withdrawn_sd` / `withdrawn_because` preserved and the paper's actual
  (qualitative) hexanal content recorded. **No substitute value invented** — none exists.
* `source_status` deliberately **NOT** set to `no_verifiable_source`: the DOI is valid, the
  paper is real and open, and the qualitative ordering *is* genuinely from it. Asserting both
  would contradict the DOI and trip citation-gate check 3.
* `source_citation_full` added with the verified citation; `transferability_notes` rewritten
  (it read *"High-value direct PPI anchor for hexanal retention"*).
* `data/lit/slr_incorporation_matrix.json` `entries[10]`: `exact_numeric_anchors` withdrawn and
  replaced by an honest statement plus the values the paper *does* support (2-octanone 14.73 →
  octanal 52.76±4.65; heptanal 13.73 → trans-2-heptenal 54.60±3.95; cis-4-heptenal 18.19 →
  trans-2-heptenal 54.60±2.98, all pooled, all from figures); `confidence_tier` **medium_high →
  low**; `incorporation_status` → `encoded_modeled_shown_pending_reanchor`; `next_action`
  rewritten to RE-ANCHOR OR RETIRE **[P, owner]**.
* Consequence worth stating plainly: this is a **wrong-compound anchor of exactly the class
  Wave K rated FATAL**, on the repo's single most-modelled analyte. It was reference-only at
  runtime (`headspace.py` reads `runtime_surrogate`, not `numeric_reference`), which bounds the
  damage but does not excuse the tier.

### (e) T1-04 / T1-05 — the Pratap-Singh origin, and the 1-hexanol factors that still ship

Wave K recorded 260 / 380 / 80 / 120 ppb as *"origin unknown, no derivation found"*. **They are
not unknown.** All four are printed in `maillard_validation_benchmarks.md` **§3.1** as
`~260 ± 35`, `~380 ± 42`, `~80`, `~120` — the tilde-hedged rows of the same table whose two
**bold** rows (2-pentylfuran `638 ± 49`, `2492 ± 199`) Wave K confirmed *verbatim*. Same file,
same commit era and same forensic fingerprint as §1.3, which Wave S2b proved fabricated. **The
fingerprint — bold/unhedged = transcribed, `~`-hedged = invented — predicts Wave K's result with
no misses**, which is its second independent validation.

* `content_correction_note` on **both** `data/benchmarks/*_40C_PratapSingh2021.json` gains
  `origin_identified_2026_08_27_wave_t3` naming §3.1 and the fingerprint. **No value changed** —
  Wave K/M/O already corrected them.
* `src/matrix_calibration_registry.py` header gains a Wave T3 block: what this changes (the
  *status* of the two 1-hexanol factors) and what it does not (nothing about the hexanal
  factors, which Wave O correctly refitted onto the verified 1138.00 / 1621.71 ppb).
* **The two 1-hexanol observability factors are now labelled at the record level.** Both the pea
  `1.0` (which *defines* the lane: `0.063 × 1269.8 = 80`) and the soy `0.143 / 0.063 =
  2.269841` get `source: "no_verifiable_source (2026-08-27, Wave T3) — back-solved from a
  fabricated value"` plus a full `notes` trail. These are **live constants**: the paper reports
  `n.d.` for hexanol in both matrices and gives soy's entire alcohol fraction as 40 ± 9 ppb of
  1-octen-3-ol — one third of the 120 ppb the numerator was solved from — and this lane carries
  the external hold-out's worst miss, **li_2026_hme 1-hexanol at 1117×**.
* **NOT REFITTED AND NOT RETIRED, and that is the whole point.** There is no measurement to fit
  to; substituting a plausible number would be the original defect again. Refit-or-retire is an
  **OWNER DECISION [P]** — see (h).

### (f) T1's blanket warnings

* **`data/benchmarks/maillard_validation_benchmarks.md` — file-level header added**, ahead of
  everything, because this file lives in the directory a reviewer trusts most and its §6
  restates its numbers as module *contracts*. It states: (1) the document is
  abstract-reconstructed, with the **bold-vs-tilde fingerprint** given as a one-sentence triage
  rule and both of its independent validations named; (2) §1.3, §3.1 and §3.2 as **confirmed**
  sources of laundered values that reached shipped code, each with what it produced and where
  that landed; (3) that §2.6, §3.2 and §4.4 carry "retrieve via WoS/Scopus" placeholder DOIs, so
  a DOI here is not evidence the paper was located; (4) that §6 inherits every defect; (5) the
  rule, verbatim from the Gemini README — **no number may enter the model from this file**; and
  (6) that the *absence* of a section-level warning box means nothing has been checked, not that
  a section is clean. §3.2's placeholder DOI is quoted with spaces inserted.
* **`docs/research/archives/README.md` — created.** Five files under a `docs/` path that reads
  as authoritative, all LLM output (2 Gemini digests, the Arrhenius "Exhaustive Analysis", 2
  Elicit reports), duplicating the `data/Gemini_Deep_Research/` corpus and carrying **no
  warning at all**. The README tables what each file actually is, restates the rule verbatim,
  adds two cautions specific to these copies (their internal citations have never been
  CrossRef-checked; digest-adjacency is not verification — with a pointer to (a)), explains why
  they are kept rather than deleted (they are the evidence trail), and records the **live
  dangling reference**: `data/species/desirable_targets.yml` and `off_flavour_targets.yml` both
  cite `pathways.md`, **which exists nowhere in this repository**, as the authority for the
  odour thresholds that are the denominator of every OAV in `src/sensory.py`. Preferred
  disposition (`git mv` under `data/Gemini_Deep_Research/`) recorded, not done.

### (g) T1-09 — the `_v1` mirror repair: 18 DOIs propagated, all CrossRef-verified

The dominant shape T1 found: a `benchmark_intake_registry.json` record that *did* get a DOI in
the 2026-08-26 repair, and its `computational_priors.json` `_v1` twin — **the copy the runtime
reads** — that did not. Scripted, diff-only, never hand-edited.

Method, deliberately conservative: candidate pairs required an **exact base-name match** after
stripping `_v1`/`_v2`/`_prior`/`_anchor`/`_point`/`_baseline`, a **single** unambiguous DOI in
the sibling's canonical `doi` field in `benchmark_intake_registry.json`, and then **every DOI
was resolved live against `api.crossref.org/works`** and its authors/year/title/journal compared
against the twin's own citation string. Propagation ran only when the citation's named surname
appears in the DOI's author list AND the years are within 1 — or when the citation names no
author at all (a bare "Ref. N" / journal-only pointer, which has nothing to contradict).

**20 candidates → 20 DOIs resolved → 18 propagated, 2 declined.** Each propagated record gets a
`doi_propagation` block recording the date, the wave, that the action is a provenance repair
with no value change, the sibling and file it came from, the full CrossRef result (title,
authors, year, container), the match basis, and why the repair was needed. The `_v1` twins now
anchored: `matoba_1988_nucleotide_hydrolysis_v1`,
`aliani_2005_donor_potency_nucleotide_context_v1`, `martins_2001_maillard_kinetics_modelling_v1`,
`van_boekel_2001_maillard_kinetics_review_v1`, `acs_jafc_0c01925_protein_binding_hierarchy_v1`,
`wang_2023_mft_retention_prior`, `mottram_2001_bmfd_retention_prior`,
`siripitakpong_2026_fft_retention_prior`, `jafc_2020_egcg_deoxyosone_trapping_v1`,
`jafc_2019_ref24_polyphenol_thiol_capping_v1`, `blank_2001_epoxydecenal_guardrail_v1`,
`frankel_1982_c182_hexanal_scission_v1`, `esterbauer_1991_4hne_kinetics_v1`,
`kamal_eldin_2003_triolein_scission_v1`, `glomb_1995_3dg_fragmentation_stoichiometry_v1`,
`frontiers_2022_hcw_aa_arrhenius_v1`, `scielo_brasil_aa_crosslink_hierarchy_v1`,
`jafc_2019_ref21_pea_gum_arabic_architecture_v1`.

**DECLINED — 2, because they are re-attributions, not propagations** (recorded rather than
performed; a wave that silently rewrites an attribution is the defect this audit exists to
remove):

1. `hidalgo_zamora_2004_4hne_pentylpyrrole_v1` cites *Hidalgo & Zamora (2004), JAFC 52:7126*;
   the sibling's DOI `10.1021/acs.jafc.5b01502` resolves to **Globisch, Kaden & Henle (2015)**.
   Neither named author appears. **[P]**
2. `grosch_1999_c183_propanal_scission_v1` cites *Grosch & Wieser (1999)*; the sibling's DOI
   `10.1007/bf02542413` resolves to **Ullrich & Grosch (1988), JAOCS 65** — 11 years and a
   different first author. **[P]**

One near-miss worth recording because it shows the rule is a rule and not a vibe:
`kamal_eldin_2003_triolein_scission_v1` was first rejected because CrossRef spells the author
`Kamal‐Eldin` with U+2010 HYPHEN and the repo uses ASCII `-`. Confirmed by Unicode
normalisation and propagated — **by normalising, not by relaxing the criterion**; the
`match_basis` field says so.

**LEFT UNANCHORED: 40 records** in `computational_priors.json` that carry no `doi` field for
their id anywhere in `data/lit/*.json` and no `source_status` flag. (T1's figure of 42 is an
*id*-level count across all of `data/lit`; 40 is the record-level count inside the one file the
runtime reads, after this wave.) **Not invented — flagged.** Each gets
`source_anchor_status: "unanchored_no_doi_field"` and a dated note stating that this is a
**disclosure, not a verdict**: it explicitly declines to assert `no_verifiable_source`, because
some of these citations may well be retrievable and nobody has tried. Retrieval is **[P]**.

### (h) T1-08 — fixed at source rather than waived (the one record check 6 caught)

`ref41_ppi_sulfur_volatile_binding_v1` cited *"DOI ref. 41 in
raw/11_maillard_lipid_crosstalk.md"* — a reference **number** inside an LLM research dump —
while claiming `provenance_tier: literature_derived_transfer`, a tier asserting a primary source
was read. The evidence for demotion was already in the tree: the only DOI ever attached to that
claim (`10.1021/acs.jafc.9b06882`, on the sibling registry record) is recorded in
`results/validation/citation_verification_ledger.md` as **TOPIC-MISMATCH** — Zhao et al. 2019, a
rodent reproductive-toxicology paper.

Tier demoted **`literature_derived_transfer` → `direct_model_assumption`** (an existing tier in
this repo's vocabulary, and the honest one: the DMTS > DMDS > DMS ordering and the binding-site
counts are assumptions this model makes, not measurements it inherited), plus
`source_status: no_verifiable_source` and a `source_status_note`. The `source` string is **kept
verbatim** rather than emptied, for two reasons: it is the evidence trail, and it is pinned as
such by `tests/scientific/test_runtime_first_registry_landing.py:121`. **No value, no
condition, no runtime behaviour changed** — verified by the Family-11 tests.

### (i) TESTS RE-PINNED — 3, each with a dated causal comment; none relaxed

1. `tests/scientific/test_slr_reference_payloads.py` — `jafc_3c05991_hexanal_binding`
   `incorporation_status` `encoded_modeled_shown` → `encoded_modeled_shown_pending_reanchor`.
   **Cause:** (d). The 16-line comment carries the verbatim quote and the reasoning.
2. `tests/scientific/test_honest_headline_guards.py::test_no_verifiable_source_census_...` —
   **102 / 80 / 62 → 120 / 98 / 80**, and the function renamed to match. **Cause:** this wave
   labelled 18 already-shipping records — 15 in `protein_source_registry.json`, 2
   `runtime_surrogate` blocks in `retention_reference_payloads.json`, 1 `ref41`. The guard's own
   docstring said counts *"must never rise silently"*; the re-pin therefore itemises which
   records moved and why, and the guard text is **strengthened**: a rise is now acceptable only
   when it is a labelling correction of numbers that were already shipping AND the re-pin
   accounts for which records moved.
3. Same file, `test_the_data_qm_records_that_moved_the_census_from_84_to_102_...` — the
   `data/lit` half of the split **84 → 102**. The `data/qm` half is **unchanged at 18**, which
   is what that test exists to protect.

**The honest reading of 62 → 80**, stated in the README, in AUDIT.md and in the test: **80 was
always the true count of unverifiable numbers the runtime consumes. 62 was an undercount.**
Nothing got worse; 18 things that were already being eaten got labelled.

### (j) DOC SYNC

* `README.md` "On literature provenance": 102 / 80 / 62 → **120 / 98 / 80**, restructured into
  the two rises (84→102 = `data/qm` un-gitignored; 102→120 = this wave's labelling), each with
  its cause, and stating explicitly that no value was added, changed or invented.
* `AUDIT.md`: the headline table row, the `data/qm` narrative paragraph, and the "The rest"
  paragraph all carry the new triple with the same account.

### (k) GATES

| Gate | Result |
|---|---|
| `scripts/ci/citation_gate.py` | **PASS** — 84 files, 991 DOI-bearing fields, 318 unique DOIs, **0 waivers**, with check 6 newly live |
| `scripts/ci/holdout_guard.py` | **PASS** — 3/3 |
| `scripts/ci/fit_target_gate.py` | **PASS** — 2/2 |

Targeted tests (no full-suite run, per the owner directive): 132 unit tests over
`test_headspace`, `test_literature_runtime`, `test_matrix_calibration_registry`,
`test_matrix_correction`, `test_matrix_recalibration`, `test_matrix_branch_deltas`,
`test_hexanal_nonanal_calibration`, `test_calibration_scope`, `test_literature_family_registry`,
`test_external_validation_report`, `test_usability_reports`; plus 82 over the scientific
provenance/matrix/audit files and 12 over the honest-headline and validation-contract guards.
All green after the three re-pins.

### (l) WHAT THIS WAVE DECLINED, AND WHY

1. **Refitting or retiring the 1-hexanol observability factors** (T1-04). No measurement exists
   to fit to; retiring the lane changes predictions and the hold-out. **Owner decision [P].**
2. **Widening `matrix_uncertainty_factor` for mock-sourced protein profiles** (T1-01). A science
   decision, not a labelling one. **[P]**
3. **The two DOI re-attributions in (g).** Attaching a DOI whose authors contradict the
   record's own citation is a content change wearing a provenance change's clothes. **[P]**
4. **T1-07 (`extrusion_r12_*`, four records citing "Literature Report 12" with no DOI
   anywhere).** Out of this wave's brief and not caught by check 6, because those records claim
   no `provenance_tier` at all. A stricter variant of check 6 — *any* record citing a digest
   must carry `no_verifiable_source` — would catch all four plus two
   `deep_research_backlog.json` entries. Implementable, deliberately not implemented here: it
   needs the six fixes to land in the same change, and this wave did not own them. **[P]**
5. **Re-sourcing the odour thresholds in `data/species/*.yml`** (T1-10). Warned about in (f)
   instead. Per-compound retrieval carries a real risk of substituting plausible numbers for
   absent ones, which is the defect the whole remediation exists to remove. **[P]**
6. **Re-deriving `log_slope` as ln(1.59)/2** (T1-02). <1% prediction change, no epistemic gain,
   breaks a behavioural pin.
7. **Deleting `scripts/trace_key_values.py` / its report.** T1 offered deletion or inversion and
   preferred inversion; the evidence is worth more than the tidiness.

### [P] CARRIED FORWARD FROM WAVE T3

1. **1-hexanol observability factors (`0.063` pea / `0.143` soy) — refit or retire.** Live
   constants back-solved from `~80` / `~120` ppb, numbers that exist only in this repo's own
   brief. Owner call; do not let a later wave quietly refit them to something plausible.
2. **`jafc_3c05991_ppi_hexanal_binding` — re-anchor or retire.** There is no PPI hexanal
   percent-bound in that paper, and its qualitative statement puts PPI at the *lowest* binding
   affinity.
3. **`protein_source_registry.json` — empirical substitution, per source, with citations**, or
   an explicit widening of `matrix_uncertainty_factor` for mock-sourced profiles.
4. **The 40 flagged-but-unanchored `computational_priors.json` records** — attempt retrieval;
   convert `unanchored_no_doi_field` into either a verified DOI or an evidenced
   `no_verifiable_source`.
5. **The 2 declined re-attributions** (`hidalgo_zamora_2004_4hne_pentylpyrrole_v1`,
   `grosch_1999_c183_propanal_scission_v1`).
6. **T1-07 / the stricter check-6 variant** (see (l)(4)).
7. **`data/species/*.yml` odour thresholds and the dangling `pathways.md` pointer** (T1-10).
8. **`git mv` `docs/research/archives/` under `data/Gemini_Deep_Research/`**, once the
   `data/species/*.yml` filename references are repaired.
9. **`results/validation/key_value_trace_report.md` is gitignored**, so the disarmed wording
   exists only on disk and in the script. If that report is meant to be reviewable, it needs a
   `.gitignore` exception like its siblings.

## Wave T4 — the ketose rearrangement escaped both pH and water-activity corrections (2026-08-27)

Executes the Tier-4 bugs from the Wave T2 round-2 dead-code survey. **One defect moves numbers
(one directional-panel row flips from MISS to OK); the rest are unreachable-anchor and
stale-prose repairs.** No correction curve was reshaped, no constant was refitted, and the
authoritative screening barrier table (`FAST_BARRIERS`) was not touched.

### (a) T4.1 — `Heyns_Rearrangement` was emitted by the engine and absent from every list that enumerates the families, including two that asserted it absent

**The defect, measured before it was touched.** `src/reaction_templates.py:60` emits
`Heyns_Rearrangement` for any ketose + amino acid; the engine produces exactly one such step
for fructose + glycine, and `tests/unit/test_smirks_engine.py::TestHeynsRearrangement::
test_heyns_fires` has asserted it fires since before Wave S1b. The family was nevertheless
missing from `ENGINE_FAMILY_LABELS` and from **all four** Wave S1b family sets, while its
aldose twin `Amadori_Rearrangement` was in two. Consequence, measured:

| pH | Amadori ion x aw | Heyns ion x aw | ratio |
|---|---|---|---|
| 5.0 | 0.000999 x 0.8286 = **0.000828** | 1.0 x 1.0 = **1.000** | **1208x** |
| 6.0 | 0.009901 x 0.8286 = 0.008204 | 1.0 | 121.9x |
| 7.0 | 0.090909 x 0.8286 = 0.075325 | 1.0 | 13.3x |

Two steps that are the same reaction on the other sugar were being corrected by a factor of
1000 differently.

**How it survived four waves — three independent blind spots, all the same shape.** Each guard
measured something narrower than the thing it reported on:

1. `tests/unit/test_wave_i_tooling.py::_engine_family_labels_in_sources` AST-scanned only
   string **literals** passed to `ElementaryStep(...)`. `reaction_templates.py:56-79` binds the
   label to a local first (`family = "Heyns_Rearrangement"`) and passes
   `reaction_family=family`. Invisible.
2. Wave S1b's family census (`_emitted_family_net_water`) enumerated `data/benchmarks/*.json`.
   **No shipped benchmark uses a ketose** — verified this wave by resolving every benchmark's
   sugar list, live *and* quarantined: **zero ketoses in either**. So the ketose branch was
   never enumerated by any guard.
3. On the strength of (2), `src/conditions.py`'s `get_ph_multiplier` docstring recorded
   *"heyns — no Heyns family is emitted"*, and
   `test_the_substring_keys_that_matched_nothing_stay_documented_as_dead` asserted it. A true
   statement about the panel, published as a statement about the engine — and the most
   dangerous line in the file, because `"heyns"` is the **only** key in branch 2 that matches
   the family (`"amadori"` is not a substring of `"heyns_rearrangement"`), so a cleanup agent
   acting on that sentence would have silently deleted the ketose rearrangement's pH dependence
   from the three ungated lanes (`kinetics.py`, `pathway_ranker.py`, `cantera_export.py`).

**What was changed.**

* `src/conditions.py` — `heyns_rearrangement` added to `_ALPHA_AMINO_NUCLEOPHILE_FAMILIES` and
  to `_LABUZA_EMPIRICAL_FAMILIES`. The aw membership is derived from **measured** net water,
  counted off the actual emitted step in the S1b idiom, not assumed:
  `OCC(=NCC(=O)O)C(O)C(O)C(O)CO -> O=CC(NCC(=O)O)C(O)C(O)C(O)CO`, one reactant, one product,
  no water either side — **net 0, uniform**, byte-identical to Amadori, so it takes no
  mass-action term and takes the empirical browning curve for the same reason Amadori does.
* `src/conditions.py:~417` — the false docstring line **corrected, not deleted**, and `"heyns"`
  relabelled LIVE AND LOAD-BEARING with the reason it must not be removed.
* `src/family_sensitivity.py` — `Heyns_Rearrangement` added to `ENGINE_FAMILY_LABELS`.

**Exactly-once discipline, asserted rather than argued.** Wiring the explicit sets does NOT
make the `"heyns"` substring redundant or double-applying, and this was checked in both lanes:
in the prediction lane `get_ph_multiplier` is reachable only through
`_enolisation_route_ph_correction`, which admits the three enolisation branch-point families
and nothing else, so Heyns takes the ionisation term once and the route term never (identical
to Amadori); in the three ungated lanes `get_ph_multiplier` is called directly and
`_ionization_correction` is never called, so the substring is the only pH term either family
gets. The import-time disjointness assertions were **extended**, not relaxed:
`_LABUZA_EMPIRICAL_FAMILIES` is now asserted disjoint from both water arms. That gap was real
and previously unasserted — `_water_activity_correction` returns from the first arm that
matches, so an overlap would not double-multiply, it would **silently shadow**, handing a
water-shedding family the peaked browning curve with no error anywhere.

**HONESTY ON THE SCIENCE.** This is a **widening**, not a de-substring-ing, and S1b's
"reproduces EXACTLY what the old substrings reached" comment no longer holds verbatim; it is
labelled as such in the source. The old `"amadori"` substring would not have matched
`heyns_rearrangement` either, so **the gap predates Wave S1b** — what S1b added was the false
claim. Further caveat, recorded in the source rather than buried: putting
`amadori_rearrangement` in the alpha-amino-nucleophile set is itself inherited from the legacy
substring, not re-derived — for a *rearrangement* the amine is already bound and the 1,2-proton
shift is rate-limiting, so the free-base fraction is a proxy for the upstream condensation
equilibrium rather than for a nucleophilic attack. **This wave makes Heyns consistent with
Amadori; it does not settle whether the shared treatment is right.** See [P] below.

**MEASURED EFFECT — stated in both directions.**

* **Benchmark panel: EXACTLY ZERO, and structurally so.** No live or quarantined benchmark
  resolves a ketose sugar, so no Heyns step is emitted anywhere in the panel. Confirmed by
  running `tests/scientific/test_benchmarks.py`, `test_primary_benchmark_campaign.py`,
  `test_wave_s1_additive_flux_2026_08.py`, `test_free_aa_quantitative_regression.py`,
  `test_validation_contract.py`, `test_honest_headline_guards.py`,
  `test_wave_p_chemistry_2026_08.py` — all green, no snapshot re-pinned. This refines the
  Wave T2 expectation ("fructose appears only in quarantined benchmarks"): it appears in
  **no benchmark's sugar list at all**, only in prose.
* **Directional panel: +1 row. `20/29 -> 21/29` agree** (8 unevaluable, unchanged);
  `sugar_identity 7/8 -> 8/8`. **SUG-12 flips MISS -> OK.**

  | row | before | after | literature |
  |---|---|---|---|
  | SUG-12 HMF, fructose vs glucose | 896.7 vs 1617.3 (**A<B, wrong way**) | **3325.4 vs 1617.3 (A>B, 2.06x)** | 1350 vs 451.5 mg/kg (2.99x) |
  | SUG-13 furfural, glucose vs fructose | 898.3 vs 842.2 (OK, 1.07x) | 898.3 vs **15.5** (OK, 57.9x) | 333 vs 276 ppb (1.21x) |
  | SUG-11 pyrazine, fructose term | 0.3318 | 0.3353 | (row unevaluable either way) |

  **Every changed number in the whole panel involves fructose; nothing else moved by a digit.**

  Mechanism, stated so the win is not over-claimed: throttling the Heyns step (~100x from
  ionisation at pH 6, ~0.3x from Labuza at aw 0.3) diverts fructose flux out of the
  Heyns -> 3-deoxyosone -> furfural branch and into `Fructofuranosyl_Dehydration`, the
  ring-retained route Wave P added on Perez Locas & Yaylayan 2008. So the *direction* of the
  improvement is the one the mechanism predicts. **But SUG-13 is now right for a more extreme
  reason than the literature supports** — the model's glucose/fructose furfural ratio went from
  1.07x to 57.9x against a measured 1.21x, i.e. that row's ordinal agreement is now bought with
  a magnitude error an ordinal panel cannot see. **Not tuned, not clawed back, and flagged**:
  it is evidence that the fructose branch partition wants a quantitative anchor, not that the
  model got better at furfural.

### (b) T4.1 blind-spot closure — the scanner AND a runtime guard

`tests/unit/test_wave_i_tooling.py`:

* `_engine_family_labels_in_sources` now resolves `reaction_family=<Name>` back to the string
  constants assigned to that name. The widening is **dataflow-scoped on purpose**: a blanket
  `family = "..."` scan produces false positives (`matrix_experiment_intake.py:179` assigns
  `family = "matrix_headspace"`, a benchmark-metadata tag, not a reaction family), so strings
  are harvested only for names that actually reach an `ElementaryStep`/`_step` call, within the
  function that does it. Pinned by a new `test_the_static_scan_now_sees_labels_bound_to_a_local_variable`,
  which asserts both that the two locally-bound labels ARE seen and that the two metadata
  strings are NOT.
* **`test_engine_family_labels_cover_runtime_emission` — the stronger form, and the real fix.**
  It parses nothing: it runs the engine over three parametrised pools and reads
  `step.reaction_family` off the emitted steps, so no scanner blind spot can hide a family. The
  **ketose pool is the point** and the test says so in its docstring: the static scan and the
  S1b census both effectively enumerated the benchmark panel, which contains zero ketoses,
  which is how this survived.
* `tests/scientific/test_wave_s1b_ph_aw_routing_2026_08.py` — `_emitted_family_net_water` now
  enumerates the panel PLUS a named `_ENGINE_ONLY_PRECURSOR_POOLS` list, so every claim built
  on that census is about the engine again. The dead-key test drops `"heyns"` from the dead
  loop and **asserts it live instead** — `== ["heyns_rearrangement"]`, i.e. strictly stronger
  than the old "no hits", because it pins which family the key reaches and that it is the sole
  matcher. Two new tests:
  `test_the_ketose_rearrangement_is_corrected_like_its_aldose_twin` (pins **equality** with
  Amadori across pH 5/6/7/9 x aw 0.3/0.8/0.95 for all three correction terms, so a future
  re-derivation must move both together) and
  `test_each_family_takes_exactly_one_water_activity_treatment` (the aw counterpart of the
  existing pH-twice guard; the shadowing failure had no assertion until now).

### (c) T4.2 — the `thiamine_degradation` Arrhenius anchor was unreachable

`data/lit/arrhenius_params.yml` has always carried a `thiamine_degradation` A/Ea pair that no
emitted family could reach: the engine's thiamine steps emit `Additive_Thermal_Degradation`,
which canonicalises to its own `FAST_BARRIERS` key, so `_arrhenius_yaml_key` returned `None`
and `cantera_export.py:157-163` silently fell back to the heuristic prefactor. One line added
to `yaml_key_map`. Verified by direct call: `_arrhenius_yaml_key("Additive_Thermal_Degradation")
== "thiamine_degradation"` and `get_arrhenius_params` returns `(1.0e9, 24.09 kcal, "dead_doi")`.
Reachability re-measured over 6 precursor pools: YAML keys no emitted family maps to went
`['dehydration', 'thiamine_degradation']` -> `['dehydration', 'mutarotation']` (`dehydration` is
documented-deliberately-inert; `mutarotation` is T4.3).

**Scope and two caveats, all recorded at the change site.** Cantera export lane ONLY;
`get_barrier("Additive_Thermal_Degradation")` is still 25.0 and no shipped prediction moves
(pinned by a new `test_wiring_the_thiamine_anchor_left_the_screening_lane_untouched`).
(1) The family is **lumped** — it covers the four thiamine steps *and* the glutathione cleavage
at `reaction_templates.py:2052`, so GSH now inherits the thiamine A/Ea; the screening lane
already lumps both at 25.0, and there is no `glutathione_cleavage` YAML entry to split them.
(2) **The anchor is unsound and this does not launder it**: `source_quality: dead_doi`, the DOI
404s, the source string hedges *"or similar lit"*, Ea = 100.8 kJ/mol is unanchored. Wiring it
makes the export report `dead_doi` where it previously reported `heuristic` for a step whose
literature entry existed all along — more informative, not more anchored. The new test asserts
`quality == "dead_doi"` precisely so a silent upgrade to `literature_estimated` fails.

### (d) T4.3 — the `Sugar_Ring_Opening` / mutarotation lane cannot fire. DOCUMENTED, NOT DELETED.

`SmirksEngine.enumerate` runs it as "Pre-Phase: Sugar Ring Opening", i.e. presents it as a live
prerequisite of the cascade. It requires a cyclic hemiacetal, and **every** sugar the precursor
resolver can produce is open-chain (and it knows only glucose / ribose / xylose / fructose).
A 30-line honest note added at `_sugar_ring_opening`, recording the measurement, the two
`FAST_BARRIERS` keys and the YAML entry that advertise the lane, and — the part worth keeping —
that the model currently **assumes** the ring/open-chain equilibrium lies fully open by shipping
open-chain SMILES, which is a modelling shortcut, not a chemical truth (D-glucose is ~99% cyclic
at equilibrium). Deleting the rule would remove the machinery and leave the shortcut
undocumented. Pinned by `test_the_mutarotation_anchor_is_still_unreachable_and_that_is_recorded`,
which fails loudly if a sugar ever resolves to a ring form. **No chemistry deleted.**

### (e) Small honesty fixes carried from Wave T2

* `_select_accumulating_projection_species` — two AST-verified unused parameters removed
  (`steps`, and `downstream_margin_kcal: float = 0.25`, a magic default a reader could have
  spent an afternoon tuning to no effect). Both call sites updated (`src/recommend.py`,
  `tests/unit/test_budget_projection.py`). Behaviour unchanged by construction.
* `debug_channel_flux` payload comment — said the inner keys are *rate-limiting `step_key`s*.
  They are the full ordered step-set (`_route_channel_id`). "Rate-limiting step_key" is the
  description of the rule that function's own docstring records as implemented, **measured and
  REJECTED**. Corrected; no behaviour change.
* `direct_sulfur_bonus` audit comment — named `Thiol_Addition_Norfuraneol` and
  `Thiol_Addition_Legacy_Shortcut` as "the routes that reach MFT today". Neither is emitted any
  more (Wave N retired the first on isotope evidence; the second has zero literals and zero
  runtime emissions). Corrected to `Thiol_Addition_Pentodiulose` /
  `Mercaptoketone_Cyclodehydration`. **The knob itself was NOT touched** — it stays deliberately
  dead per its own documentation, and the conclusion is now *stronger* than when Wave H wrote it.

### (f) Curated-layer drift — dated PARITY NOTE added, chemistry NOT added

`src/curated_pathways.py` gains a header note recording, from an independent count made this
wave: **Wave N mirrored** (`Deoxyosone_Reduction` + `Thiol_Addition_Pentodiulose` present, the
retired norfuraneol step gone); **Wave P NOT mirrored** — none of the six Wave P families appear
(the C2+C3 mercaptoketone lane to MFT, `Fructofuranosyl_Dehydration`,
`Furanone_Reductive_Opening`, `Furanone_Amino_Acid_Reduction`); and the ketose lane absent too.
The Wave P ledger section contains zero occurrences of "curated", so the mirror was not skipped
by decision — it appears not to have been considered. **Consequence recorded in the file:
`scripts/generate_reaction_network.py` draws `docs/assets/reaction_network.pdf` from this layer,
so the published architecture figure depicts pre-Wave-P chemistry.** No chemistry was added:
hand-adding six families to a mirror nothing tests is how the drift got here. Owner call below.

### (g) OPEN ITEMS ADDED

- [P] **Re-derive the aldose/ketose rearrangement pH treatment, or confirm it.** Wave T4 made
  `Heyns_Rearrangement` consistent with `Amadori_Rearrangement`; it did not establish that the
  shared treatment is correct. For a rearrangement the amine is already bound, so the
  Henderson-Hasselbalch free-base fraction is a proxy for the upstream condensation
  equilibrium, not for a nucleophilic attack. Both families must move together — pinned by
  `test_the_ketose_rearrangement_is_corrected_like_its_aldose_twin`. Calibration work, owner
  sign-off.
- [P] **The fructose branch partition wants a quantitative anchor.** After T4.1 the model's
  glucose/fructose furfural ratio is 57.9x against a measured 1.21x (SUG-13), while its HMF
  ratio 2.06x sits reasonably against 2.99x (SUG-12). The ordinal panel scores both as
  agreements. Not tuned this wave.
- [P] **`Sugar_Ring_Opening` / mutarotation (T4.3)**: implement cyclic sugar forms in
  `data/species/` + the resolver — at which point the rule, `FAST_BARRIERS["mutarotation"]`,
  `["ring_opening"]` and the YAML `mutarotation` entry all become live and the open-chain
  assumption becomes explicit chemistry — OR retire the lane and all four together. Do not do
  the latter piecemeal.
- [P] **Curated layer (f)**: mirror Wave P and regenerate `docs/assets/reaction_network.pdf`, OR
  declare the layer deliberately frozen at the Wave N topology and say so in the figure caption.
  The durable fix is probably a **parity test** between the two family vocabularies rather than
  another hand-sync.
- [P] **Two incompatible path-span formulas in `src/recommend.py`. FLAGGED, NOT UNIFIED.**
  Wave S1 extracted the step arithmetic into `_evaluate_step`, which computes the path span as
  the log-sum-exp series resistance (`src/recommend.py:~1173`). Two blocks below, the same
  method recomputes a path barrier by hand with the **old max rule**:
  `src/recommend.py:1526` (lipid trapping efficiency) and `src/recommend.py:1568` (lysine / DHA
  budget), both `path_barrier = max(max_r_dist, barrier)`. They feed
  `metrics["trapping_efficiency"]` and `metrics["lysine_budget_dha"]`, **not** `predicted_ppb`,
  so no headline number is affected — but the module carries two definitions of "path span",
  one of which this ledger calls the physically motivated one. Consolidating changes two
  published PBMA metrics; ~40 lines of duplication. **Deliberately not unified blind.**
- [P] **Re-anchor `thiamine_degradation` (T4.2)** against the real thiamine kinetics the intake
  registry already names (Voelker 2018 `10.1016/j.foodres.2018.06.056`; Voelker 2021
  `10.1186/s13065-021-00773-y`), and consider a separate `glutathione_cleavage` YAML entry so
  the lumped `Additive_Thermal_Degradation` family stops giving GSH the thiamine pair.

### (h) GATES

| Gate | Result |
|---|---|
| `scripts/ci/citation_gate.py` | **PASS** — 84 files, 991 DOI-bearing fields, 318 unique DOIs, 0 waivers |
| `scripts/ci/holdout_guard.py` | **PASS** — 3/3 |
| `scripts/ci/fit_target_gate.py` | **PASS** — 2/2 |

Targeted tests only, per the owner directive (**no full-suite run**): **268 passed, 0 failed**
across `test_wave_i_tooling`, `test_wave_s1b_ph_aw_routing_2026_08`, `test_arrhenius_params`,
`test_budget_projection`, `test_conditions`, `test_wave_h_2026_08`, `test_smirks_engine`,
`test_chemistry_soundness`, `test_data_integrity`, `test_uncertainty_propagation`,
`test_wave_s1_additive_flux_2026_08`, `test_free_aa_quantitative_regression`,
`test_family_sensitivity`, `test_validation_contract`, `test_honest_headline_guards`,
`test_wave_p_chemistry_2026_08`, `test_benchmarks`, `test_primary_benchmark_campaign`,
`test_pentose_hexose_sulfur_ordering`, `test_pipeline`. **11 tests added this wave.**
The uncommitted work of Waves S1b, S2c and T3 was preserved throughout; nothing was committed
or stashed. Orchestrator certifies the batch.

## Wave U — the reaction network's FIRST out-of-sample test (2026-08-27)

Wave S1's structural finding is the reason this wave exists: **all four bundles in the
external hold-out run the `matrix_only` execution path.** They read a lipid-oxidation load off
a matrix profile and return before `Recommender.predict_from_steps` is ever called. That is
why Waves S1, S1b and S2c each moved dozens of in-panel rows and left **all eight hold-out
points bit-identical** — three consecutive waves, one of them a shipped *barrier* revert. The
invariance was evidence about the hold-out's coverage, not about the model. **The chemistry
this repository is about had never been scored on a system it had not already seen.**

This wave harvested twelve content-verified free-precursor literature points that DO exercise
the network, froze them under the hold-out guard, and recorded the CURRENT model's predictions
so the pending rate-calibration wave can be scored against predictions committed **before it
ran**. The ordering is the deliverable.

**Nothing was tuned to these points, by this wave or anyone.** No barrier, projection
constant, observability factor or prior was touched. `src/` was not edited at all except
`scripts/ci/holdout_guard.py`, which was *strengthened*.

### (a) THE POINTS — 12 bundles, 22 targets, every value quoted from a retrieved source

All live in `data/benchmarks/external_validation/maillard_path/`, all
`evidence_class: external_validation_only`, all `maillard_path_holdout: true`, all
`execution_path: free_precursor`, all `protein_type: free`.

| # | System | Source (access) | Target(s), native unit | Dedupe |
|---|---|---|---|---|
| 1-4 | **ribose 25 mM + cysteine 25 mM**, 0.5 M K-phosphate pH 5.5, sealed 20 mL Duran, 3 mL; 100 °C/4 h, 110 °C/2 h, 120 °C/1 h, 130 °C/0.5 h | Yiltirak, Balagiannis, Koek, Koch & Elmore 2026, *Food Res. Int.* 231:118600, `10.1016/j.foodres.2026.118600` — **full text + publisher supplementary** (CentAUR PDF via `--resolve` to 46.22.140.67; Table S3 unzipped out of Elsevier `mmc1.docx`) | MFT 6.88 / 3.29 / 2.4 / 1.71 µg/L; FFT 1.28 / 1.46 / 1.68 / 1.62 µg/L — **SIDA** (MFT-d3, FFT-d2, matrix-matched curve) | **0 hits** for the DOI or either lead author across `data/ docs/ src/ results/` |
| 5-6 | **glucose 0.25 M + L-alanine 0.25 M**, unbuffered, sealed ampoule, 130 ± 1 °C, 2 h, at **pH 5** and **pH 8** | Schibilsky 2019, TU Berlin dissertation (DepositOnce bitstream `969d7075-…`) — **full text PDF**, Tabelle 10 + Tabelle 17 read directly | furfural 17/25 µM; furaneol 9/46 µM; HMF 454/798 µM — internal standard with **measured** response factors (Tabelle 10) | **0 hits**; no DOI exists, typed `identifier`/`identifier_scheme` pair used (same convention as the Liu 2021 NC State thesis) |
| 7 | **glucose 0.2 M + asparagine 0.2 M**, 0.1 M phosphate pH 6.86, 4 mL sealed tube, 180 °C, 30 min | Ye et al. 2024, *Foods* 13(17):2836, `10.3390/foods13172836` — **full text** (Europe PMC `PMC11395303` fullTextXML) | acrylamide **140.58 ± 13.92 µmol/mol Asn** — 13C3-acrylamide **isotope dilution** LC-QQQ-MS/MS. Kept in the paper's own basis-free ratio unit; **no ppb conversion** | **0 hits** |
| 8-10 | **glucose 27.75 mM + asparagine 33.3 mM**, pH 6.0, 100 mL, 180 °C; 10 min and 30 min in 1% acetic acid, plus 30 min in deionised water | Chang / Sung et al. 2021, *Polymers* 13(12):1901, `10.3390/polym13121901` — **full text** (`PMC8229482`) | acrylamide 28 / 1459 / 832 ppb; HMF 7 ppm — authentic-standard calibration curves | **0 hits** |
| 11 | **fructose 27.75 mM + asparagine 33.3 mM**, pH 6.0, 100 mL, 180 °C, 30 min | Lin, Ting, Ndraha, Hsiao & Sung 2022, *Polymers* 14(8):1565, `10.3390/polym14081565` — **full text** (`PMC9031937`) | acrylamide 1859 ppb; HMF 12.28 ppm | **0 hits** |
| 12 | **glucose 555 mM alone** (10% w/v), no amino acid, pH 4.36 measured, 121 °C, 18 min autoclave | Steinhagen et al. 2021, *Pharmaceuticals* 14(11):1121, `10.3390/ph14111121` — **full text** (`PMC8625795`), Table 4 read directly | 5-HMF 17.4 ± 3.9 µg/mL — LC-MS/MS, calibration 0.5-100 µg/mL, R² ≥ 0.993 | **0 hits** |

**Diversity, since it was the stated criterion:** sulfur thiols (8 targets), furanone +
oxygen-heterocycle (6), acrylamide (5), HMF (4); temperatures 100-180 °C; times 10 min - 4 h;
pH 4.36 - 8.0; three different sugars and four different amino acids; one **temperature/time
series**, one **fixed-temperature time course**, one **pH pair**, one **sugar-identity pair**
across bundles 10/11, and one **sugar-only** system with no amino acid at all.

**One methodological note, recorded because it is a choice.** Bundles 7-11 set
`metadata.family` to something other than `"safety"`. `_evaluate_loaded_benchmark` dispatches
`family == "safety"` to the lumped `predict_acrylamide` solver, which never touches the
network. Routing around that dispatch is the entire point of the bundle, and each says so in
its `metadata.notes`; the safety solver's answer for the same system is recorded alongside as
`cross_path_reference` and is explicitly not the scored number.

### (b) DECLINED CANDIDATES — every one, with the reason

Nine sources were retrieved far enough to judge and then refused. Two of the refusals are the
most useful output of this wave.

| Candidate | Retrieved to | Reason declined |
|---|---|---|
| **Martins 2003**, WUR thesis, glucose/glycine 0.2 M pH 6.8, **HMF max 20 µmol/L at 120 °C** | full text (`edepot.wur.nl/121418`) | **CONTAMINATED — this is a calibration source.** `data/lit/arrhenius_params.yml` takes Ea = 97 kJ/mol from that chapter's Table 5.2, and `data/lit/timeseries/martins2005_glucose_glycine_80_100_120C_pH68.yml` holds **ten species curves digitised off its figures**. 10 grep hits. Scoring the same chapter's HMF as a "hold-out" would have been laundering. This is why the independent Schibilsky thesis was used for the furanone lane instead. |
| **Knol 2005 / Knol 2008 thesis**, `10.1021/jf050504m`, acrylamide 1.5 mol% (glucose) and 3 mol% (fructose) | full text (`edepot.wur.nl/122035`) | Three reasons, any one sufficient: the yields are stated in concluding prose with **no temperature or time attached** (they are maxima over a 120-200 °C × 0-45 min series); the per-timepoint concentrations exist **only as figures**; and Table 7.2's `a` column is a **fitted logistic asymptote**, not a measurement. Also partially contaminated — `FAST_BARRIERS["safety_risk_acrylamide"] = 30.83` is a declared transfer of the Knol group's Ea. Used in this wave **only** to source bundle 7's reactant molarity, which is recorded as an inference. |
| **Blank, Fay, Lakner & Schlosser 1997**, `10.1021/jf960997i`, HDMF 2.6-5.1 µg / EHMF 6.8-10 µg per mmol pentose | abstract only (ACS 403; EPFL Infoscience 429 on every attempt) | The yields are a **RANGE across unspecified conditions**, and the precursor concentrations are not in the abstract. Taking an interior point of a published band is precisely the 342/200 failure mode Wave S2b exposed. 3 existing grep hits for the DOI. |
| **Hofmann & Schieberle 1998**, `10.1021/jf9705983`, MFT 1.4 mol% | abstract only (closed on Unpaywall, OpenAlex, S2, Europe PMC, IA Scholar) | **The mol% denominator is defined nowhere in the abstract**, and neither are the precursor concentrations or the reaction volume. This is the paper Waves S2b/S2c retired; re-ingesting its abstract number would have re-created the retired anchor under a new name. |
| **Grosch, Zeiler-Hilgart, Cerny & Guth 1993**, MFT 0.2 µg/L (cys+ribose) vs 300 µg/L (cys+thiamine) | second-hand, quoted verbatim in Cerny 2015 (chapter retrieved in full) | Second-hand, and the conditions are the single word **"boiling"** — no temperature, time, pH, buffer or vessel. Usable as a ~1500× ordinal claim about the thiamine route, not as an executable point. **Worth an ILL: it is the only absolute thiamine-route number found anywhere.** |
| **Cerny & Briffod 2007**, `10.1021/jf062874w` | abstract only | Closed everywhere. **Correction to this wave's own brief:** the assumption that it "reports quantified µg amounts of MFT/FFT as a function of pH" **could not be confirmed** — the published abstract contains only 13C-label distributions and relative statements, and no absolute unit. `data/lit/slr_incorporation_matrix.json` already carries a `citation_caveat` saying the repo's numbers do not come from it; that caveat stands and is now independently corroborated. |
| **Ping et al. 2021**, `10.3390/antiox10070993`, acrylamide **2250 ± 82 µg/kg**, true 13C3 SIDA | full text (`PMC8300766`) | **The most painful decline.** Genuine isotope dilution, clean control, verified table. But the system is 10 µmol each of asparagine and glucose **homogenised into 1 g of silica gel**: the paper states neither a water activity nor a reaction volume, so executing it requires inventing **two** free parameters, and its µg/kg basis refers to a solid support rather than the aqueous medium the projection layer assumes. Measured at HEAD 12f43dd, the prediction swings **47× across a plausible aw range** (70.1 ppb at aw 0.2 → 1.5 ppb at aw 0.9). Declined rather than scored with two invented parameters. *The one robust statement it does support, recorded here and not in a bundle: the model is at least **32× low** on this system at every aw tested.* |
| **Chen et al. 2026**, `10.1016/j.crfs.2026.101433`, five pyrazines in µg/L with per-compound calibration curves (R² > 0.995) | full text (`PMC13186052`) | Not a free-precursor system: the precursor is a **purified Thr-Amadori compound**, and the engine has no ARP entry point. Substituting threonine + xylose would be modelling a different experiment. |
| **Ratsch et al. 2025**, `10.3390/molecules30030535`, SIDA-quantified furfural / 2-methylpyrazine with deuterated standards, fructose 0.02 M + one amino acid 0.01 M | full text (`PMC11821130`) | The vehicle is **Chardonnay base wine** — an undefined precursor background plus ethanol. Closest thing in open access to an absolute furanone-family model system, and still not one. |

Rejected outright on quantification method, each verified by reading the paper's own sentence:
**Wang et al. 2021** `10.3390/foods10020273` (peak-area normalisation, RF assumed 1, and *"m
represents the sample mass (g)"* never defined — a ~50× basis ambiguity; 5 existing grep hits);
**Zhang et al. 2022** `10.1016/j.ultsonch.2022.105913` (the running text calls the same numbers
both "peak area" and "concentration", and 2,3-DMP ≫ 2,5-DMP by ~1000× is chemically
implausible); **Cha et al. 2019** `10.1038/s41598-019-41824-8` (α-dicarbonyls *are* absolute
but none is in `data/species/`; the volatiles are *"peak area ratio (PAR)"*); **Song et al.
2019** `10.3390/polym11030521` (*"assuming that the relative response factor was 1 and the
recovery ratio was 100%"*); **Ma et al. 2022** `10.3389/fnut.2022.940202` (SIDA method, but
**no absolute concentration anywhere in text or tables** — figures only; DOI already in repo);
**Scalone 2016** Ghent thesis (*"GC-MS peak Area × 10⁶"* for every model system; its only
SIDA-grade assay is on cookies and bread); **Zheng 2023** `PMC9868338` and **Foods 2023**
`PMC9818664` and **Okumura 1990** (single IS with RF = 1, or "semiquantified ... NIST library",
or "Area %").

**ILL asks this wave generated** (all CrossRef-verified, all closed on every open route
exercised — Europe PMC, Unpaywall, OpenAlex, Semantic Scholar, IA Scholar, CentAUR,
DepositOnce, mediaTUM, biblio.ugent.be, J-Stage, EPFL Infoscience):
`10.1021/jf9705983` **Tables 1-6** (the mol% precursor matrix — would settle the retired
anchor); `10.1021/jf0200826` Mottram & Nobrega; `10.1021/jf980980v` + Whitfield & Mottram 2001
(the H₂S + norfuraneol branch); `10.1021/jf062874w` Cerny & Briffod (to settle whether it holds
absolute numbers at all); `10.1021/jf990954c` Hofmann, Münch & Schieberle 2000 (**the** Strecker
paper — nothing open substitutes for it); Grosch et al. 1993 in *Progress in Flavor Precursor
Studies*. A separate finding worth carrying: across ~15 Europe PMC open-access query
formulations, **no open-access paper reports absolutely-quantified Strecker aldehydes from a
free-precursor sugar + amino-acid model system.** The Strecker branch cannot be held out from
open sources.

### (c) THE FROZEN ARTIFACT — `results/validation/maillard_path_holdout_frozen_predictions.{json,md}`

Generated by `scripts/generators/generate_maillard_path_holdout_frozen_predictions.py` at
**git HEAD 12f43dd** (branch `audit-remediation`), before any rate-calibration wave saw these
points. The file records that commit in its own header and states that S3 must be scored
against **it**, not against a fresh run. **It was un-gitignored in the same change**
(`.gitignore`, with the reason written in place): an untracked pre-registration is not a
pre-registration, because nothing would stop a later wave regenerating it after calibrating
and comparing it with itself.

| Measure | Result |
| --- | --- |
| Bundles / targets | 12 / 22 |
| Quantitatively scored | 21 |
| **Median fold error** | **6.0388×** |
| Median abs log10 error | 0.7809 dex |
| Worst / best | **52.586×** (FFT, 130 °C) / **1.5168×** (MFT, 100 °C) |
| Within 10× | **12 / 21** (57.1%) |
| Within-point orderings | **8 / 12** pairs |
| Cross-bundle response directions | **3 / 6** |
| Structural zeroes | **1** |

**Unit handling.** Only two non-ppb target units are scored quantitatively, and both are
basis-free ratios: `mol_percent` and `umol_per_mol_limiting_precursor`. The predicted value is
derived from the model's OWN declared identity and nothing else —
`ppb = molar × MW × ppb_conversion_factor (=1e6)` from `src/projection.py` — divided by the
lowest declared precursor molarity, which is the same limiting-precursor rule
`_estimate_projection_budget` uses. Bundle 7 additionally states that the paper's denominator
(mol of asparagine charged) and the model's denominator coincide, which is what makes the pair
commensurable at all. The µM → ppb conversions in bundles 5-6 use only a molecular weight and
the thesis's own volumetric basis, and the exact arithmetic is written into each target's
`native_unit_and_conversion` field. **No conversion in this wave introduces a molar basis.**

**`structural_zero` is counted separately and deliberately.** Bundle 12 (glucose alone) gets a
prediction of exactly 0.0 against a measured 17 400 ppb. It has no finite fold error, so it
cannot enter the median — which is precisely why it is named and counted rather than dropped
as "not evaluable". A median that silently omits total misses reports a better model than
exists.

### (d) WHAT THE BASELINE ACTUALLY SAYS — the median is the least interesting number

**The encouraging part is real and should be said first.** Wave S2c ended with *the sulfur
branch now has zero absolute literature anchors*. This wave gave that branch its first real
absolute measurement, as a hold-out, and it predicts **MFT at 100 °C to 1.5168×**. An
uncalibrated branch landing inside a factor of two on an isotope-dilution number it had never
seen is not nothing.

**Then the shape, which is wrong in three specific ways a fold-error median cannot show:**

| Series | Axis | Compound | Measured | Predicted | Response ratio | Direction |
|---|---|---|---|---|---|---|
| Yiltirak 2026 | 100 °C/4 h → 130 °C/0.5 h | MFT | ×0.249 | ×4.55 | **18.3** | **WRONG WAY** |
| Yiltirak 2026 | same | FFT | ×1.27 | ×5.48 | 4.33 | OK |
| Chang 2021 | 10 min → 30 min at 180 °C | acrylamide | ×52.1 | ×1.24 | **0.0238** | OK |
| Schibilsky 2019 | pH 5 → pH 8 at 130 °C | furfural | ×1.47 | ×0.774 | 0.527 | **WRONG WAY** |
| Schibilsky 2019 | same | furaneol (DMHF) | ×5.11 | ×2.48 | 0.485 | OK |
| Schibilsky 2019 | same | HMF | ×1.76 | ×0.774 | 0.441 | **WRONG WAY** |

1. **The sulfur branch has the temperature dependence backwards.** Measured MFT *falls* 4.0×
   across the series; the model has it *rising* 4.55×. The excellent 1.52× at 100 °C and the
   12.1× miss at 130 °C are the same error seen from opposite ends — which means a
   calibration that only shifts the sulfur barriers up or down will improve one and worsen the
   other. **This is the sharpest single instruction this wave has for S3.**
2. **Acrylamide barely responds to time.** ×52.1 measured against ×1.24 predicted, a response
   ratio of 0.024. And note *what this is not contaminated by*: the Knol Ea transfer sets the
   temperature dependence, and this series holds temperature fixed.
3. **Two of three pH responses point the wrong way**, and both wrong ones (furfural, HMF) move
   by exactly the same factor ×0.774 — i.e. they are being driven by one shared pH term, not
   by two independent ones. Furaneol, which moves correctly, is the exception.
4. **The MFT > FFT ordering fails at all four sulfur points (0/4).** All eight of the wave's
   within-point ordering hits come from the acrylamide/HMF and oxygen-heterocycle bundles,
   which score 8/8.

**A fifth finding, found while here and not looked for.** The two acrylamide lanes in this
repository disagree by **~480×** on the same 0.2 M glucose/asparagine system: the network
returns 1142.6 ppb and the lumped `predict_acrylamide` returns 544 870.3 ppb. Both are shipped.
Nothing in the panel had ever evaluated them side by side, because every acrylamide benchmark
in `data/benchmarks/` carries `family: "safety"` and therefore only ever exercises one of them.
Recorded in bundle 7's `cross_path_reference`. **[P] for the owner — not fixed here.**

### (e) GUARD COVERAGE — extended, and negative-tested

The bundles live in a **subdirectory**, `external_validation/maillard_path/`, and the reason is
scientific rather than cosmetic: `get_holdout_benchmark_files()` globs `*.json`
**non-recursively**, so these points do not enter `external_validation_report`'s median. That
median is a matrix-only lipid-oxidation number, and averaging two different execution paths
into one figure is the kind of thing this campaign exists to stop. It also means no
`test_honest_headline_guards` pin had to move.

That placement would have left them **entirely unguarded**, so `scripts/ci/holdout_guard.py`
was extended — strengthened only, nothing relaxed:

* **Check 1 now uses `rglob`, not `glob`.** It reports 16 bundles (4 legacy + 12 new); before
  the change it would have seen 4 and said PASS.
* **New check 4, `check_maillard_path_bundles`.** Every file flagged `maillard_path_holdout:
  true` must (a) live under the hold-out directory, (b) carry
  `evidence_class: external_validation_only`, (c) declare
  `metadata.execution_path: free_precursor` — a bundle that silently drifted to `matrix_only`
  would keep scoring and stop testing anything — and (d) have its `benchmark_id` absent from
  the whole fit/calibration corpus (`data/lit/**/*.json` plus every
  `results/validation/*refit*.json` / `*rederivation*.json`).

**Negative-tested, three ways, each confirmed to exit 1 with the check-4 line suppressed:**
planting `{"fit_target": "mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026"}` in
`results/validation/`; flipping one bundle's `execution_path` to `matrix_only`; and copying a
flagged bundle up into `data/benchmarks/`. All three were reverted; the clean tree exits 0.

`tests/scientific/test_wave_u_maillard_path_holdout.py` adds **38 assertions** covering the
same invariants from the test side, plus one that every target carries a `target_quote`, a
`target_value` and a stated basis, and one that the frozen artifact still names its git
commit. **The tests deliberately do NOT pin the scores** — pinning the number you just
produced is the circularity Rounds 1-3 removed. The one pinned number is the bundle count, so
that dropping an inconvenient point breaks the build.

### (f) [P] CARRIED FORWARD

- [P] **Score S3 against `maillard_path_holdout_frozen_predictions.json` at HEAD 12f43dd.**
  Not against a regenerated copy. The file says so itself.
- [P] **The two acrylamide lanes disagree by ~480×** (see (d)). Owner decision: reconcile them,
  or state in the README which one is authoritative for which regime.
- [P] **The sulfur branch's temperature dependence has the wrong sign.** A uniform barrier
  shift cannot fix it (see (d)(1)). Whatever S3 does here, the 100 °C and 130 °C points must
  be reported together.
- [P] **Two pH responses move by the identical factor ×0.774**, suggesting furfural and HMF
  share a single pH term where the measurement says they should not. Worth a look at
  `_enolisation_route_ph_correction` before calibrating anything.
- [P] **The free-precursor network returns exactly zero for a sugar with no amino-acid
  partner.** Either wire a caramelization lane or state the scope limit on the surfaces that
  advertise a `carbohydrate_pyrolysis_and_caramelization` family.
- [P] **The Strecker branch cannot be held out from open sources** (see (b)). It needs
  `10.1021/jf990954c` through an institutional route or it stays untested.

### (g) GATES

| Gate | Result |
|---|---|
| `scripts/ci/holdout_guard.py` | **PASS** — **4/4** (was 3/3), 16 bundles checked, 12 flagged, 33 fit records scanned |
| `scripts/ci/citation_gate.py` | **PASS** — 96 files (was 84), 1001 DOI-bearing fields, 323 unique DOIs, 0 waivers |
| `scripts/ci/fit_target_gate.py` | **PASS** — 2/2 |

Targeted tests only, per the owner directive (**no full-suite run**):
`tests/scientific/test_wave_u_maillard_path_holdout.py` **38 passed**. No existing test was
re-pinned, because no shipped number moved: this wave added data, one generator, one gitignore
negation and one guard check, and edited no module under `src/`.

---

## Wave S3 — the first RATE-LEVEL calibration of the sugar trunk (2026-08-27)

Every barrier in this repository had, until this wave, been either an endpoint fit, a
literature midpoint, or an estimate by analogy. Nothing had ever been fitted to a measured
**trajectory**. Wave Q's 767 point-verified concentration-time values made that possible for
the first time; this wave spends them.

### (0) PRE-REGISTERED PREDICTION — written before the hold-out was re-scored

The mission required predicting which of Wave U's three named structural errors a trunk rate
calibration could move, *before* measuring. Recorded here first, at the point in the wave
where the fit existed and the hold-out had not been re-run. Scored in section (e).

The predictions follow from one mechanical fact established by reading the screening lane
end to end (`src/barrier_constants.py:804` → `benchmark_validation.py:662-673` →
`recommend.py:1179, 1211-1224, 910, 952` → `projection.py:170-223`):

> **The FAST screening lane consumes barriers ONLY as branching ratios.** The predicted
> total ppb is `total_volatile_budget_molar × mole_fraction × MW × 1e6`, and
> `total_volatile_budget_molar` comes from `projection.py:170-223`, which never sees a
> barrier — it is built from the limiting precursor molarity, the duration, and a *separate*
> apparent Ea of 120 kJ/mol. The flux is additionally normalised twice
> (`recommend.py:910` by `max_weight`, `:952` by `total_activity`). Measured: applying a
> uniform +2 / +5 / +10 kcal/mol shift to every barrier changes the total predicted ppb by
> under 0.4%, all of it molecular-weight re-weighting.

| # | Wave U structural error | Prediction | Confidence |
|---|---|---|---|
| 1 | **Sulfur T-dependence backwards** (MFT measured ×0.249 across 100→130 °C, model ×4.55) | **Not fixed.** The temperature dependence of the magnitude is set by `projection.py`'s own 120 kJ/mol apparent Ea, which no trunk parameter touches. Trunk barriers can shift the sulfur branching ratio versus its competitors, and that ratio is T-dependent, so the number may *move* — but nothing in a sugar-trunk fit knows the sign of the sulfur error. | ~85% not fixed; direction of any movement not predicted |
| 2 | **Acrylamide ~40× under-responsive in time** (×52.1 measured, ×1.24 predicted over 10→30 min) | **Cannot be fixed, structurally.** The budget is `1 − exp(−k_app·t)` with `k_app·t ≈ 3e-3 ≪ 1`, i.e. **linear in t**. The ceiling on the model's response to a 3× time change is therefore ×3, against a measured ×52. No barrier value anywhere can lift a linear response to a 52-fold one. | ~95% |
| 3 | **Furfural and HMF share one pH term** (both move by exactly ×0.774) | **Unchanged.** The entire fit corpus sits at a single pH (6.8). No pH parameter is fitted, and the shared term lives in `_ionization_correction` / `_enolisation_route_ph_correction`, which this wave does not touch. | ~95% exactly unchanged |
| — | **Median fold error (6.0388×)** | Moves little. Magnitudes are budget-driven and barrier-invariant to first order; only branching ratios can move, and only where a route crosses the `_temporal_accessibility` saturation threshold. | predicted \|Δ median\| < 20% |
| — | **The structural zero** (glucose alone → 0.0 vs 17 400 ppb) | **Stays zero.** No amino acid means no condensation step, so the trunk never starts. A rate calibration cannot create a caramelization lane that does not exist. | ~99% |
| — | **Directional panel** (21/29; pH 4/7, aw 0/3) | pH and aw rows exactly unchanged (no pH or aw parameter fitted). Headline moves by at most ±1. | ~90% |

**The honest summary of that table: this wave was predicted, in advance, to be a calibration
that the hold-out cannot see.** That prediction is itself the finding — it says the
repository's absolute accuracy problem does not live in its barriers.

### (a) LANE — the fit does not live in either shipped lane, and that was measured

Three candidate lanes existed. Each was tested by running it, not by reading it.

**The FAST screening lane cannot be fitted to a trajectory because it has no time axis.**
Tracing `src/barrier_constants.py:804` → `benchmark_validation.py:662-673` →
`recommend.py:1179, 1211-1224, 910, 952` → `projection.py:170-223`: there is no ODE solver
anywhere in `src/` (`grep -rn "solve_ivp\|odeint\|scipy.integrate" src/` returns nothing).
Two findings fell out of that trace and both matter later:

* **The per-family prefactor cancels exactly.** `conditions.get_rate_constant` multiplies by
  `A(family)`; `effective_barrier_from_rate_constant` on the next line divides by the same
  `A(family)`. Measured: `bar_eff = Ea − RT·ln(multiplier)`, identically, for every family.
  Every `A_value` in `data/lit/arrhenius_params.yml` is **dead in the screening lane**.
* **The one surviving absolute rate is family-agnostic.** `recommend.py:1216` calls
  `get_reference_pre_exponential()` with **no argument**, so every family gets the 1e11
  default regardless of its YAML entry.

**The Cantera export lane integrates ODEs but cannot represent the experiment.** Four
blockers, all measured at HEAD b2a1b20 on glucose + glycine:

1. **THE INITIAL MOLARITY IS DISCARDED.** `src/kinetics.py:373` does
   `phase.X = initial_concentrations`; Cantera normalises `X` to sum to one, so only the
   *ratio* survives and the absolute scale comes from Girolami molar-volume estimates.
   **Measured: feeding 0.02 M, 0.2 M and 2.0 M gives bit-identical trajectories** (Amadori
   at 9000 s = 1.93585 kmol/m³ in all three), against a t=0 glucose concentration of
   2.546 mol/L that nobody asked for. A bimolecular constant cannot be fitted through that.
2. **Four of the ten measured species do not exist in the network** — no fructose (no Lobry
   de Bruyn isomerisation step at all), no formic acid, no acetic acid, no melanoidins.
3. **The YAML lumps the enolisations.** `_arrhenius_yaml_key` maps `1,2-enolisation`,
   `2,3-enolisation` and `enolisation_intermediate` onto one `enolisation` A/Ea pair — three
   separately-measured steps forced to share one rate.
4. **The YAML barrier is not even the operative number.** `cantera_export.py:161` does
   `barrier_kcal = max(barrier_kcal, float(barrier_lit))` — a max against the FAST barrier,
   not an assignment — after which every step is marked reversible with its equilibrium set
   by Joback-estimated NASA polynomials.

**So the fit lives in a third, dedicated lane: `src/trunk_kinetics.py`.** Explicit
mass-action ODEs, `solve_ivp`, native units (mmol/L, minutes, kJ/mol), no thermo estimation,
no gating. Eight steps, six species. Distorting a shipped lane to make it fittable would
have been worse than the disease. **No module under `src/` imports it**, and
`test_wave_s3_trunk_kinetics.py` asserts that by scanning the tree, so the day it becomes
false a test says so.

### (b) FITTED PARAMETERS — every one with its CI and its constraint status

Corpus: `martins2005_glucose_glycine_80_100_120C_pH68.yml` (glucose, DFG, 3-DG, 1-DG at
80/100/120 °C) plus the **pH 6.8 series of `martins2003_DFG_amadori_degradation.yml`** — an
experiment on *isolated* DFG that separates the Amadori compound's decomposition from its own
formation, which is what makes the branch split identifiable at all. **176 values, 16 free
parameters.** Multi-response weighted least squares on `log(y + floor)`; floors from each
file's own stated read-off uncertainty; **sigmas estimated from the replicate scatter in the
data itself** (glucose 0.016, DFG 0.040, 3-DG 0.076, 1-DG 0.118 in log space), not assumed.
Arrhenius carried as (k at 100 °C, Ea) — Martins' own reparameterisation, because over an
80-120 °C window A and Ea are too correlated for a fitted A to mean anything.

| step | repo family | k(100 °C) | 95% CI | Ea kJ/mol | 95% CI | status |
|---|---|---:|---|---:|---|---|
| `k_schiff` | `schiff_condensation` | 1.71e-05 L/(mmol·min) | 1.48e-05 – 1.97e-05 | 94.4 | 86 – 103 | **well constrained, both** |
| `k_amadori` | `amadori_rearrangement` | 0.1536 /min | 0.104 – 0.226 | 34.7 | 13 – 57 | constrained |
| `k_dfg_3dg` | `enolisation_intermediate` | 6.46e-03 /min | 3.96e-03 – 1.05e-02 | 144.8 | 116 – 174 | constrained |
| `k_dfg_1dg` | `2,3-enolisation` | 3.05e-02 /min | 5.58e-03 – 0.167 | 73.0 | 15 – 131 | weak (30× band) |
| `k_dfg_other` | *(none — structural gap)* | 4.31e-03 /min | 2.95e-08 – 630 | 89.5 | **−314 – 493** | **SLOPPY / Ea sign unconstrained** |
| `k_3dg_out` | `1,2-enolisation` (lump) | 0.1047 /min | 0.0568 – 0.193 | 135.6 | 96 – 176 | constrained |
| `k_1dg_out` | `furanone_amino_acid_reduction` (lump) | 3.031 /min | 0.534 – 17.2 | 40.1 | **−21 – 101** | weak / Ea sign unconstrained |
| `k_glc_other` | *(none — structural gap)* | 2.71e-04 /min | 4.31e-05 – 1.70e-03 | 146.1 | 26 – 266 | weak |

**Three honesty items about that table, none of them cosmetic.**

1. **Reduced χ² = 14.69.** The sigmas are the data's own replicate scatter, so this says the
   model's systematic error runs **3.8× the measurement noise**. Every CI above has been
   widened by that factor; unscaled Gauss-Newton intervals would have been 3.8× too tight
   and would have called five of these parameters "well constrained".
2. **An Ea confidence interval that reaches below zero is not a wide interval, it is no
   interval.** A negative Arrhenius Ea is unphysical, so `k_dfg_other` and `k_1dg_out` have
   *no determined temperature dependence*, and the report says exactly that rather than
   printing a tidy number. The three sloppiest Hessian eigendirections (eigenvalues 2.3e-05,
   2.8e-04, 1.3e-03 against a condition number of 5.6e+08) are **all** Ea directions —
   `Ea k_dfg_other`, `Ea k_glc_other`, and the `Ea k_1dg_out`/`Ea k_dfg_1dg` pair. The rates
   are what this data constrains; several of the activation energies are not.
3. **Only 1 of 60 random starts reached the optimum, and the first version of this fit
   reported the wrong one.** The Schiff/Amadori profile re-optimises 15 of 16 parameters at
   each of 20 fixed ratios, which sweeps the surface far more widely than the multi-start —
   and it found a point **better than the reported optimum by Δcost 7.2**, i.e. the profile
   would have been published with negative Δcost entries in it. Rather than tolerance it
   away, the generator now adopts any such point, re-polishes all 16 parameters from it and
   recomputes the profile, terminating only when a full profile finds nothing better
   (`profile_informed_repolish_rounds: 1` in the artifact). This objective is rugged and the
   record says so.

**Fit quality, per response group:** median fold error 1.01–1.24× on all twelve
glucose/glycine groups and 1.24–1.38× on the two isolated-DFG groups.

**Two steps correspond to NO repository reaction family**, and are declared as gaps rather
than folded into a neighbour. `k_glc_other` lumps the glucose→fructose isomerisation and the
amine-independent fragmentation to organic acids — the repository's SMIRKS network emits
**neither**, and at 100 °C they account for a large share of the measured glucose loss
(glucose falls 70 mmol/L while glycine falls only 29). `k_dfg_other` is the DFG channel that
does *not* regenerate glycine.

### (c) CROSS-CHECK — one independent, one that only looks independent

The distinction is load-bearing and is enforced in the artifact's own text.

**Brands (2002), thesis Chapter 4 — GENUINELY INDEPENDENT.** Different amine (protein-bound
lysine ε-amino, not free glycine), different sugar:amine ratio (10:1, not 1:1), different
author, different data set, recovered by Wave S2. None of it entered the objective.

| step | fitted @120 °C | Brands @120 °C | ratio | fitted Ea | Brands Ea |
|---|---:|---:|---:|---:|---:|
| condensation (sugar + amine → Amadori) | 8.04e-05 L/(mmol·min) | 2.4e-04 ± 2e-05 | **0.33×** | 94 | 114 ± 2 |
| **TOTAL Amadori degradation** | 0.189 /min | 0.2805 | **0.67×** | — | 126 |

**This is the wave's credibility test and it passes.** Two laboratories, two amines, two
independent fits, agreeing to **1.5× on the Amadori degradation rate** and 3× on the
condensation. The degradation row is the stronger of the two, because the fitted side is
anchored by an isolated-DFG experiment Brands never ran.

**Martins (2003) M4 — NOT INDEPENDENT, and labelled so everywhere it appears.** Martins
fitted his constants to the very data fitted here; agreement is a *reproducibility* result —
a different implementation, a different scheme and a different objective recovering the same
rates — and is not evidence from a second experiment.

| fitted step | Martins step(s) | fitted k(100 °C) | comparator | ratio |
|---|---|---:|---:|---:|
| `k_schiff` | 1 | 1.71e-05 | 1.61e-05 | **1.06×** |
| `k_dfg_3dg` | 4 | 6.46e-03 | 1.11e-02 | 0.58× |
| `k_dfg_1dg` | 7 | 3.05e-02 | 1.57e-02 | 1.94× |
| `k_dfg_other` | 6 | 4.31e-03 | 7.08e-03 | 0.61× |
| `k_3dg_out` | 5 + 9 | 0.1047 | 0.1969 | 0.53× |
| `k_1dg_out` | 8 | 3.031 | 1.45 | 2.09× |
| `k_glc_other` | 2 | 2.71e-04 | 1.64e-03 | 0.17× |

**A units bug in this very table, caught and fixed, recorded because the repository has been
burned by exactly this class of error.** Two of Martins' M4 steps are **bimolecular** — step
1 (Glu + Gly → DFG) and step 9 (3-DG + Gly → melanoidins) — so their X values are in
L/(mmol·min). The first draft compared step 9's 8.12e-4 directly against a first-order
constant and reported a spurious **131× disagreement**; converting it to pseudo-first-order
at the experiment's own 200 mmol/L glycine gives 0.53×. The generator now converts a
bimolecular comparator **only** when the fitted step it faces is first-order, and records
that it did.

### (d) THE SCHIFF/AMADORI AUTHORITY INVERSION — RESOLVED. NEITHER FILE WAS RIGHT.

Open as `[P]` since Wave I: `FAST_BARRIERS` says the Schiff condensation is ~2.0e4× faster
than the Amadori rearrangement; `arrhenius_params.yml` says the Amadori rearrangement is
~3.3e4× faster than the Schiff condensation — a disagreement of ~6.6e8 in the ratio about
which of the first two steps of the entire cascade is rate-determining.

**The Martins data measure the Amadori intermediate directly, and the isolated-DFG
experiment measures its decomposition alone. The data decided.**

The two steps can only be compared as pseudo-first-order rates, because the condensation is
bimolecular and the rearrangement is not. At the experiment's own 200 mmol/L glycine, 100 °C:

```
    k_schiff × [Gly]  =  3.42e-03 /min      (condensation)
    k_amadori         =  1.54e-01 /min      (rearrangement)
    ratio             =  44.9               (Amadori faster)
```

A profile likelihood over that ratio — with `k_amadori` pinned and all 15 remaining
parameters re-optimised at each point, from four starts each — gives a **95% interval of
40–45**. Every decade from 0.1 to 1e6 outside it is rejected:

| claim | ratio asserted | verdict | χ² statistic |
|---|---:|---|---:|
| `FAST_BARRIERS` (screening lane) | 5.0e-05 | **WRONG SIGN** | ~7.5e+04 |
| `arrhenius_params.yml` (Cantera lane) | 3.3e+04 | right sign, **wrong by ~700×** | ~9.1e+02 |
| **fitted** | **44.9** (40–45) | — | 0 |

**So the reconciliation is not "pick a file".** The YAML has the ordering right and the
magnitude wrong by about three orders of magnitude; the screening lane has the ordering
backwards. **The condensation is rate-determining and the rearrangement is fast** — which is
also why Martins was entitled to lump them into one step, and why the fitted `k_schiff`
reproduces his lumped step 1 to 1.06×.

Written into **both** files' headers, with the arithmetic. The YAML header additionally
carries the repair constraint its `amadori` entry has been missing: any fixed A/Ea pair must
satisfy `A·exp(−Ea/(R·373.15)) = 2.52e-03 1/s`, and with the entry's current Ea = 59.0 kJ/mol
that implies **A = 4.6e5 1/s, not the 1.0e11 written there** — i.e. the *prefactor* is the
broken half, by 5.3 orders of magnitude, exactly as that entry's own `audit_flag` suspected.
Neither file's values were changed.

### (e) SCORING THE PRE-REGISTRATION — the prediction in (0), now measured

Scored against `results/validation/maillard_path_holdout_frozen_predictions.json` at commit
`12f43dd`, **read, never regenerated**. New artifacts (append-only, the frozen file
untouched): `results/validation/maillard_path_holdout_S3_prepost.{json,md}`.

**AS SHIPPED: 0 of 22 targets moved. Every prediction is bit-identical to the
pre-registration. The directional panel is byte-identical too (21/29; `diff` clean).**

| # | prediction | outcome |
|---|---|---|
| 1 | sulfur T-dependence not fixed | **CORRECT** — MFT response ratio 4.551 unchanged (measured 0.2485) |
| 2 | acrylamide time response cannot be fixed | **CORRECT** — 1.24 unchanged (measured 52.11) |
| 3 | shared pH term exactly unchanged | **CORRECT** — furfural and HMF both 0.7743 |
| — | median fold error moves <20% | **CORRECT** — 6.0388× → 6.0388×, exactly |
| — | structural zero stays zero | **CORRECT** |
| — | panel pH/aw unchanged | **CORRECT** — 4/7 and 0/3, identical |

### (f) PROPAGATION — derived in full, deliberately NOT applied, and the counterfactual measured

A principled mapping exists. The screening lane's own inverse is
`Ea = −0.001987·T·ln(k/A)`, and `A` must be the family-agnostic **1e11**, because that is
what `recommend.py:1216` uses for every family. At T = 373.15 K:

| step | family | derived kcal/mol | shipped kcal/mol | Δ |
|---|---|---:|---:|---:|
| `k_schiff` | `schiff_condensation` | **excluded** | 15.00 | — |
| `k_amadori` | `amadori_rearrangement` | **23.20** | **23.00** | **+0.20** |
| `k_dfg_3dg` | `enolisation_intermediate` | 25.55 | 21.00 | +4.55 |
| `k_dfg_1dg` | `2,3-enolisation` | 24.40 | 28.00 | −3.60 |
| `k_3dg_out` | `1,2-enolisation` | 23.49 | 28.00 | −4.51 |
| `k_1dg_out` | `furanone_amino_acid_reduction` | 20.99 | 28.00 | −7.01 |

Worked example: `k_amadori` = 0.1536 /min = 2.560e-03 /s;
`Ea = −0.001987 × 373.15 × ln(2.560e-03 / 1e11) = 23.20 kcal/mol`.

**`schiff_condensation` is excluded outright**: it is bimolecular, in L/(mmol·min), and there
is no dimensionally valid conversion to a first-order barrier without inventing a
standard-state concentration. Inventing one is precisely the 342/200 basis failure Wave S2b
exposed.

**NOT APPLIED. `FAST_BARRIERS` is unchanged by this wave.** Two reasons, the second measured:

1. **The dominant uncertainty is an assumed prefactor, not the fitted rate.** Setting
   A = 1e11 s⁻¹ for every step folds the true prefactor, the activation entropy and general
   acid/base catalysis by the phosphate buffer into Ea at **1.71 kcal/mol per decade** at
   100 °C — wider than every confidence interval in section (b).
2. **Applying it makes the model worse on both out-of-sample surfaces.** Measured, with the
   derived barriers patched in memory and never written to `src/`:

| surface | as shipped | counterfactual |
|---|---|---|
| hold-out median fold error | **6.0388×** | **7.6110×** (26% worse) |
| hold-out targets moved | 0/22 | 21/22 — 11 better, **10 worse** |
| hold-out series directions | 3/6 | 3/6 |
| hold-out within-point orderings | **8/12** | **6/12** |
| hold-out within 10x | 12/21 | 13/21 |
| directional panel headline | **21/29** | **18/29** |
| panel pH / sugar_identity | 4/7 · 8/8 | **2/7** · **7/8** |

The counterfactual is not uniformly bad — `within_10x` rises 12 to 13 and the best single
point improves from 1.52x to 1.12x — but the median, the within-point orderings and the
directional panel all move the wrong way. **On the two surfaces that were designed to be
out of sample, the verdict is consistent: do not apply.**

**One correction to my own reasoning in (0), owned rather than buried.** I argued the
magnitudes are "barrier-invariant to first order" and predicted a median move under 20%. That
is true of the *total* predicted ppb — a uniform barrier shift moves it under 0.4% — but the
per-compound ppb are branching ratios, and the derived set is not a uniform shift: it
compresses the enolisation family from a 7 kcal/mol spread to 4.6. The counterfactual moved
the median 26%. **The shipped prediction (0 movement) was right; the mechanism I gave for it
was too strong.**

**One unlooked-for confirmation of a Wave U hypothesis.** Wave U inferred, from furfural and
HMF both moving by exactly ×0.774, that they "are being driven by one shared pH term, not by
two independent ones". Under the counterfactual both move to **exactly 0.3008** — still
identical to each other after a large and *differential* barrier change. That is direct
positive evidence the shared term is a single multiplier, not a coincidence.

### (g) OTHER HONEST RESIDUALS — measured, out of the objective, not fitted to

* **Glycine yield.** The isolated-DFG experiment's Table 4.1.1 measures 83.4% / 78.1% of DFG
  returning its glycine at 100 / 120 °C. The fit — which never saw this response — predicts
  **89.5% / 90.1%**. The trunk over-attributes DFG decay to the two glycine-releasing
  channels by 6-12 points, i.e. `k_dfg_other` is fitted too small. That parameter is also the
  sloppiest in the fit, and these two facts are the same fact.
* **The 1-DG branch flux is not credible against the acid data.** Cumulative
  DFG → 1-DG flux over 150 min at 100 °C is **45.4 mmol/L** against a measured acetic acid of
  **13.6 mmol/L** — 3.3× over budget. This is exactly why `k_dfg_1dg` carries a 30× CI: the
  1-DG pool level constrains only the *ratio* `k_dfg_1dg / k_1dg_out`, and nothing in the fit
  corpus pins the pair. **Fitting the acid responses would fix it, and the repository has no
  formic or acetic acid species to fit them with.**
* **Glucose loss over-predicted 15%** (80.3 vs ~70 mmol/L at 150 min, 100 °C).
* **The missing pH physics, quantified.** The pH 5.5 series was held out of the objective
  because the model has no pH term. Median fold error against it: **6.5× at 100 °C and 14.5×
  at 120 °C.** That is the size of the gap Wave S2's 2/7 pH score was pointing at, now with a
  number on it.
* **The read-off floor, measured from the repository's own data.** The 100 °C experiment
  appears in two figures of the same thesis, extracted independently. Their median
  disagreement is 0.51 mmol/L for glucose and 0.145 for DFG — a floor no fit can beat.

### (h) WHAT CHANGED, AND WHAT DID NOT

**Changed:** `src/trunk_kinetics.py` (new, imported by nothing shipped);
`scripts/generators/generate_trunk_rate_calibration.py` and
`generate_trunk_calibration_holdout_prepost.py` (new);
`results/validation/trunk_rate_calibration_refit.{json,md}` and
`maillard_path_holdout_S3_prepost.{json,md}` (new);
`tests/scientific/test_wave_s3_trunk_kinetics.py` (new, 29 assertions);
headers of `src/barrier_constants.py` and `data/lit/arrhenius_params.yml` (the inversion
resolution — comments only, **no value moved**); `data/lit/timeseries/README.md` (its
"nothing here is wired into the calibration" claim is now false and says so);
`scripts/ci/fit_target_gate.py` + `src/fit_target_index.py` (strip `.yml`/`.yaml` as well as
`.json` from a declared fit target, kept in lockstep).

**NOT changed:** any `FAST_BARRIERS` value, any `arrhenius_params.yml` value, any projection
constant, observability factor, prior or benchmark. **No shipped number moved**, which the
0/22 hold-out and byte-identical panel confirm independently of the claim.

**The fit corpus is now machine-declared and machine-guarded.** The artifact carries
`fit_target_files` and `fit_leverage` (16 parameters / 176 rows → `global_low_leverage`), and
`fit_target_index.fit_target_map()` now resolves both YAML corpora to real ids. The `_refit`
suffix in the filename is **load-bearing, not cosmetic**: `holdout_guard` check 4 scans
`results/validation/*.json` filtered on the substring `refit`/`rederivation`, so the artifact
named `trunk_rate_calibration.json` as originally specified would have been invisible to both
gates. It now scans 34 fit records where Wave U scanned 33.

### (i) [P] CARRIED FORWARD

- [P] **Apply the derived barriers, or keep the two-lane split.** The full table is in (f).
  The measured counterfactual says applying them costs 26% on the hold-out median and 3 rows
  on the directional panel, so the recommendation is **do not apply** — but the derivation is
  on the record and the decision is the owner's.
- [P] **`arrhenius_params.yml`'s `amadori` prefactor is wrong by 5.3 orders of magnitude**,
  and the fitted rate now says what a repair must reproduce (see (d)). Do not patch Ea alone.
- [P] **The repository's network has no glucose→fructose isomerisation and no organic-acid
  species.** Both are needed before the trunk's branch fluxes can be pinned (see (g)); the
  1-DG channel is currently 3.3× over the acid budget and nothing in-corpus can catch it.
- [P] **`amadori_rearrangement` = 23.0 kcal/mol is now independently corroborated** (derived
  23.20 from a measured trajectory). This is the first FAST_BARRIERS value with trajectory
  evidence behind it, and the only one.
- [P] Wave U's `[P]` items are untouched by this wave: the ~480× acrylamide-lane
  disagreement, the sulfur temperature sign, the shared pH term (now *confirmed* shared, see
  (f)), the caramelization structural zero, and the closed Strecker branch.

### (j) GATES

| Gate | Result |
|---|---|
| `scripts/ci/holdout_guard.py` | **PASS** — 4/4; 16 bundles, 12 flagged, **34** fit records scanned (was 33) |
| `scripts/ci/citation_gate.py` | **PASS** — 97 files, 1002 DOI-bearing fields, 323 unique DOIs, 0 waivers |
| `scripts/ci/fit_target_gate.py` | **PASS** — 2/2, with the new YAML corpus declared and resolving |

Targeted tests only, per the owner directive (**no full-suite run**).

### (k) REGENERATION SCOPE — nothing downstream needed it, and that is checkable

The mission asked which artifacts consume the changed parameters. **None, because no shipped
parameter changed.** The only edits to `src/` are (i) a new module nothing imports and (ii)
comment blocks in `src/barrier_constants.py` and a two-line suffix list in
`src/fit_target_index.py`. That claim is not asserted, it is measured twice, independently:
the hold-out is **0/22 moved, bit-identical**, and the directional panel is **byte-identical
(`diff` clean)** against a run captured before any edit. Both would have caught a value that
moved. So `benchmark_summary`, `prediction_uncertainty`, `external_validation_report`, the
gap heatmap and the family figures were all left alone deliberately, and re-running them
would have produced churn with no content.

Four artifacts were **un-gitignored** (`.gitignore`, with the reason written in place).
`trunk_rate_calibration_refit.{json,md}` has to be tracked for a reason beyond review: it is
the machine-readable fit-corpus declaration that `fit_target_gate`, `fit_target_index` and
`holdout_guard` check 4 all read. Untracked, CI would run with this wave's fit **undeclared
and unguarded** — precisely the failure mode the fit-target gate exists to prevent. The
pre/post artifact is tracked because a wave whose headline is "the hold-out did not move"
owes the reader the file that shows it did not, target by target.

---

## Wave S4 — measurable binding physics against the fitted observability fudge (2026-08-27)

Wave S3 proved the absolute-accuracy deficit does not live in the barriers. This wave went
one layer down the matrix lane, to the OBSERVABILITY factors that convert a predicted total
into the number a paper reports. Every one of them was fitted; two of them
(`0.143 / 0.063`, the 1-hexanol pair) were back-solved from values Wave T3 proved appear in
no publication. **Nothing in this wave is fitted. No shipped number moved. The hold-out
entered nothing.**

### (a) THE HEADLINE, STATED BEFORE THE PLEASANT PART

**The fitted observability factors score WORSE on the never-fitted hold-out than applying no
observability factor at all**, and they beat everything in-panel. That is the overfitting
signature the mission predicted, measured:

| mode | hold-out median fold | median \|log10\| | CI coverage | in-panel median fold |
|---|---:|---:|---:|---:|
| `fitted_factors` (SHIPPED, unchanged) | **93.684x** | 1.9717 | 3/8 | **1.0004x** |
| `unit_observability` (null model, factor 1.0) | **67.419x** | 1.8288 | 4/8 | 5.9196x |
| `binding_physics` (measured constants, 0 fitted) | 68.179x | 1.8337 | **5/8** | 5.9196x |

Monte Carlo n=200 seed 0, the shipped scoring path. The `fitted_factors` row reproduces the
shipped headline **exactly** (93.6837x, 3/8 in `external_validation_report.json`), which is
the independent check that the mode machinery is inert by default.

**AND THE BINDING PHYSICS IS A WASH.** The null model is why that can be said rather than
guessed. Against `unit_observability` the binding mode moves the median by **1.1%** — because
on **12 of 16** hold-out rows there is no usable binding datum and the mode reduces to the
null model *exactly*. On the four rows where it genuinely applies (all Liu 2023, the only
water-calibrated lane in the hold-out) the two scored ones split 1-1:

* `liu_2023_ppi / nonanal` **94.22x -> 1.852x** — the repository's sharpest lipid-lane
  over-prediction against a directly-quantified reference, very nearly fixed by a constant
  measured in 2023 by someone else.
* `liu_2023_ppi / hexanal` **48.24x -> 120.9x** (against the null; 11.17x under the fitted
  mode) — worse, and for a reason that is the same reason the Pratap-vs-Liu test below fails.

Without the null model, "binding_physics 68x beats fitted 94x" would have been reported as
evidence for the binding physics. It is mostly evidence against the fitted factors.

### (b) THE SHARPEST RESULT IS A UNIT ARGUMENT, NOT A SCORE

An observability factor is the fraction of a total that the measurement sees. **It cannot
exceed 1.** The shipped constants are **4.31725 (pea ambient hexanal)** and **9.54007 (soy
ambient hexanal)**, with 5.92 on the soy furan class anchor.

So whatever those numbers are, they are not observability. They are absorbing an
absolute-scale deficit that lives UPSTREAM, in `MATRIX_BENCHMARK_BASE_MARKER_YIELDS` — and
those yields are themselves `260 / 638 / 80 / 150` divided by a common scale, i.e. built from
the very values Waves K/M/T3 retired. This is why the binding mode loses 4.3x on the
in-panel hexanal rows the moment the fitted factor is removed: the budget under-predicts
verified Pratap-Singh hexanal by exactly that factor, and the "observability factor" was the
only thing hiding it. **No binding model can repair an absolute-scale error from the
observability side.** Boundary, as the mission required: this wave did NOT touch the
volatile-budget / allocation machinery (`src/projection.py`, the base marker yields,
`src/lipid_oxidation.py`). That is the separate [P] workstream, and this result is a
quantified argument for it.

### (c) THE MODEL — one equation, zero fitted parameters

`src/protein_binding.py` (new, 0 shipped-path effect by default):

```
    K_eff  = a_p * Pow                    [L per gram of protein]
    f_free = 1 / (1 + K_eff * c_p)        [single-site Langmuir, dilute-ligand limit]
```

`a_p` is the Harrison & Hills / Viry hydrophobic interaction parameter as fitted to plant
proteins by Snel et al. 2023, and it is the ONLY per-gram binding parameter this campaign
found: every other retrieved constant is per MOLE of protein and none of those sources states
the protein molar mass it used, so the conversion cannot be done without inventing one. Those
are recorded under `not_usable_without_protein_molar_mass` rather than converted.

**The second half of the physics is the quantification basis, and it matters as much as
f_free.** `f_free` is only the right factor if the reference reports the FREE concentration:

* **matrix-matched** (standards spiked into the slurry) measures the TOTAL -> factor **1.0**,
  and `f_free` must NOT be applied. Pratap-Singh 2021 is this case, verbatim: *"The pea
  protein sample was spiked with 1, 2, 3, 4, 5 uL of stock standard and 5 uL of internal
  standard hexanal-d12 was used to generate the standard curve for quantification"*.
* **water-calibrated** (external curve built in water) reads back the FREE concentration ->
  factor `f_free`. Liu 2021 is this case, verbatim: *"pea protein solutions were replaced
  with 5 mL DI water and were spiked with known compounds at five different levels"*.
* **unknown** -> nothing is applied and the row says so. Six of the eight hold-out rows are
  this, which is most of why the binding mode has so little to do.

**NOT MODELLED, deliberately:** no denaturation term (see (f)), no covalent term, no
low-water-activity lane (an aqueous-phase Langmuir does not describe a dry extrudate).

### (d) THE ZERO-PARAMETER CROSS-CHECK — the model reproduces measurements it was not built from

`a_p` was fitted by one laboratory on 2-alkanones in demineralised water by APCI-MS. These
rows are other laboratories, other methods, and in two cases another chemical class. Nothing
was adjusted:

| record | compound | c_p g/L | measured % bound | predicted | residual (pts) |
|---|---|---:|---:|---:|---:|
| Heng 2005, pea vicilin | 2-heptanone | 10.0 | 19.00 | 19.27 | **+0.27** |
| Wang 2015, pea isolate pH 7 | 2-octanone | 8.27 | 31.90 | 32.64 | **+0.74** |
| Wang 2015, pea isolate pH 7 | 2-heptanone | 8.27 | 13.90 | 16.49 | +2.59 |
| Wang 2015, pea isolate pH 7 | 2-hexanone | 8.27 | 7.15 | 4.72 | -2.43 |
| Barallat-Perez 2023, plant isolates | octanal | 10.0 | 52.76 | 55.61 | +2.85 |
| Heng 2005, pea vicilin | 2-octanone | 10.0 | 33.00 | 36.95 | +3.95 |
| Heng 2005, pea vicilin | 2-hexanone | 10.0 | 13.00 | 5.66 | -7.34 |
| Wang 2015, pea isolate pH 8 | 2-octanone | 8.27 | 23.76 | 32.64 | +8.88 |
| Wang 2015, pea isolate pH 8 | octanal | 8.27 | 61.87 | 50.88 | **-10.99** |
| *Barallat-Perez 2023 (range endpoint)* | 2-octanone | 10.0 | 14.73 | 36.95 | *+22.22* |
| *Barallat-Perez 2023 (range endpoint)* | heptanal | 10.0 | 13.73 | 33.28 | *+19.55* |

**Median \|residual\| 3.95 percentage points over all 11; 3.27 over the 8 that are not pooled
range endpoints.** The two italicised rows are the LOW endpoints of "an increase from X to Y
for PPI, SPI and LPI", i.e. the least-binding of three proteins, not a mean — they are the two
largest residuals and they are reported separately rather than dropped.

**Two things this cross-check earns that no score could.** First, `a_p` was fitted on ketones
and both octanal rows land within 11 points, so the ketone->aldehyde class transfer is
defensible (it is nonetheless declared an UPPER bound on `f_free`, because aldehydes also bind
covalently). Second, the compound-transfer variable is validated independently: the aldehyde
series `K(nonanal)/K(hexanal)` is **28.6x** (Barallat-Perez 2025 Klotz) and **26.2x** (Liu
2026 inverse LC) against a Pow ratio of **33.1x**; and O'Keefe's 1988 glycinin
`K(1-hexanol)/K(hexanal)` is **470/270 = 1.74** against a Pow ratio of **1.78** — a 2%
agreement from a 1988 experiment that knew nothing about this model.

### (e) THE PRATAP-vs-LIU TEST — a clean, falsifiable, NEGATIVE result

Wave R's insight was that Liu's DI-water-calibrated values are lower bounds "by roughly the
bound fraction". That is now a number, and it does not do what was hoped.

```
  f_free(hexanal, pea, c_p = 100 g/L)      = 0.3990
  predicted water-calibrated under-read    = 2.5065x
  Pratap-Singh (matrix-matched, total)     = 1138 ug/L
  Liu Table 2.7 band                       = 2445 - 52454 ug/L, geometric mid 11 325
  gap BEFORE the binding correction        = 9.95x
  gap AFTER  the binding correction        = 24.94x
```

The direction is right — a water-calibrated measurement must read low — but the correction
moves Liu's numbers UP, i.e. **it widens the disagreement it was supposed to explain, by
2.5x.** The verdict is therefore that Pratap-Singh and Liu differ by MATERIAL AND LOT, not by
method: Liu's own band spans 21.45x across nine commercial lots, which swamps a 2.5x method
term. This also explains the one hold-out row where the binding physics makes things worse
(`liu / hexanal`): the model already under-predicts there, and `f_free < 1` can only push it
further under.

### (f) DENATURATION — the sources disagree in sign, so there is no term

The mission allowed an adjustment only if a source quantifies it. Three do or nearly do, and
they do not agree:

| source | material | direction | quantified? |
|---|---|---|---|
| Liu 2026 (inverse LC, 85 C / 30 min) | commercial SPI on silica | **decreases** | yes: hexanal **~-10%** in K, series -6 to -18% |
| Heng 2005 | purified pea vicilin | **decreases** | yes: octanal 96% -> 32% bound |
| Barallat-Perez 2023 | commercial isolates | **increases** | no, qualitative only |

The best-matched source (commercial SPI) gives a ~10% change in K, which moves `f_free` by
under one percentage point at every loading in the file — smaller than the Pow convention gap
alone (see (g)). A term that is disputed in sign and, where measured on the right material,
negligible, is not a term. All three are recorded in
`data/lit/binding_constants.yml -> denaturation_effect_evidence` with `modelled: false`, and
`src/protein_binding.describe_model()` reports `denaturation_modelled: false` with the reason.

### (g) THE DATA FILE — `data/lit/binding_constants.yml`, 23 records / 11 sources

Every record carries native units, the measurement conditions, a **verbatim quote**, the
table/figure location, a `verification_status`, and (where a DOI exists) a CrossRef check.
A test asserts all of that. Sources retrieved and read in full this wave:

| source | what it gave | usable? |
|---|---|---|
| **Snel et al. 2023**, Heliyon 9:e16503 (PMC10245154) | `a_p` in **L per g protein**: yellow pea **25e-5**, soy **16e-5** (ketones, R2 = 1.00); esters 11e-5 / 4.8e-5 | **YES — the only per-gram parameter found anywhere** |
| **Barallat-Perez et al. 2023**, JAFC 71:20274 (PMC10739987) | octanal **52.76 +/- 4.65** %bound and three more, at 1 wv% / pH 7 | yes, cross-check |
| **Wang 2015**, U. Manitoba thesis (mspace) | Table 5.1 pea-isolate ketone %bound across pH 3-11; Table 7.2 octanal **61.87 +/- 0.46** | yes, cross-check |
| **Heng 2005**, WUR thesis 10.18174/121674 | pea vicilin %bound, heated and non-heated | yes, cross-check |
| **Barallat-Perez et al. 2025**, npj Sci Food 9:174 | SPI Klotz K and n for hexanal (3.5e3 M^-1), nonanal (1e5), the whole C6-C10 aldehyde and 2-alkanone series | **no** — M^-1, no protein molar mass stated |
| **Liu et al. 2026**, Food Chem 525:150473 (edepot CC-BY) | SPI aldehyde K_B by inverse LC, hexanal 2623 M^-1; **the denaturation number** | **no** for K; yes for (f) |
| **O'Keefe 1988**, Iowa State dissertation | purified glycinin / beta-conglycinin K for hexanal, octanal, **1-hexanol**, 2-nonanone at 3 temperatures | **no** — M^-1 and the molar-mass sentence is ambiguous |
| PubChem / Hansch 1995 | experimental log P for hexanal (1.78), 1-hexanol (2.03), 2-octanone (2.37), 2-heptanone, 2-hexanone; XLogP3 for nonanal, 2-pentylfuran, octanal, heptanal | yes, with the estimated ones flagged |

**The dominant stated uncertainty is the Pow convention, and it is quantified.** `a_p` was
fitted against KOWWIN v1.68 estimates; the experimental Hansch values are used where they
exist. On the one compound where both are available (2-octanone: KOWWIN 2.22, Hansch 2.37)
they differ by **0.15 dex = 1.41x in Pow**, which propagates linearly into `K_eff`. For
nonanal and 2-pentylfuran no experimental value is indexed at all and XLogP3 is used, flagged
`computed_xlogp3` — a 0.3 dex error there is a 2x error in `K_eff`, and nonanal's f_free is
0.0197, so that row is the most Pow-sensitive in the file.

### (h) WIRING — a mode, not a replacement, with the no-double-count guard NEGATIVE-TESTED

Four modes, resolved by an explicit context manager, then the
`MAILLARD_MATRIX_OBSERVABILITY_MODE` environment variable, then the shipped default:
`fitted_factors` (default, unchanged), `unit_observability` (null), `binding_physics`, and
`binding_physics_out_of_domain`.

**`binding_physics_out_of_domain` is currently a NO-OP and is not scored, and the reason is
worth recording:** every non-aqueous lane in this repository (Bi roasted pea, Li 2026 HME)
also has no declared protein loading and no established quantification basis, so the mode has
nothing to compute with. **What blocks it is missing data, not the domain rule.**

THE NO-DOUBLE-COUNT GUARD. In binding mode `src/headspace.py` returns BEFORE
`get_matrix_calibration_record` is consulted, and it also drops `dynamic_release_factor` —
the second omission is the one that is easy to miss, because `compose_dynamic_retention`
routes through `resolve_compound_matrix_retention`, whose `volatile_retention` is documented
as *"fraction escaping matrix (rest is bound)"*, i.e. it is ITSELF an unanchored binding
model. Running both would count protein binding twice.

`src/benchmark_validation.py` then asserts it. The first version of the assertion demanded
`net == f_free x pH` point-by-point and **fired immediately** — correctly: the Monte-Carlo
propagator legitimately wraps that method with a sampled multiplier. The shipped form uses
the invariant that survives the sampler and still catches a leak: the sampler's multiplier is
GLOBAL to a sample while every factor being guarded against is PER COMPOUND, so the ratio
`net / (f_free x pH)` must be **identical across all four compounds on a lane**.
Negative-tested both ways in `tests/scientific/test_wave_s4_protein_binding.py`: a global
0.25x leak passes and is reported as `binding_no_double_count_ratio = 0.25`; a per-compound
3x leak on hexanal alone raises `AssertionError: ... Double counting refused`.

### (i) THE FABRICATED-FACTOR ROWS

`0.143 / 0.063` is 120 ppb over 80 ppb, neither of which appears in Pratap-Singh et al.

| benchmark | mode | observability factor | predicted ppb | measured | fold |
|---|---|---:|---:|---:|---:|
| `pea_isolate_40C_PratapSingh2021` / 1-Hexanol | fitted | 1.0 | 80.101 | *(retired — paper says n.d.)* | — |
| `soy_isolate_40C_PratapSingh2021` / 1-Hexanol | fitted | **2.2698** | 119.92 | *(retired)* | — |
| `soy_isolate_40C_PratapSingh2021` / 1-Hexanol | binding | 1.0 | 52.831 | *(retired)* | — |
| `li_2026_hme` / 1-hexanol | fitted | **2.2698** | 22 390 | 20.04 | **1117x** |
| `li_2026_hme` / 1-hexanol | binding | 1.0 | 9 866 | 20.04 | **492x** |

Deleting the fabricated factor halves the log error on the only scored row it touches, and
the two in-panel rows it also drives are unscorable because Wave K removed their target (the
paper reports n.d.). **Neither this wave nor any before it has found a pea or soy 1-hexanol
percent-bound or per-gram constant.** The nearest published numbers are O'Keefe's purified
soy fractions (glycinin 470 M^-1, beta-conglycinin 181 M^-1 at 293 K) and Wongprasert 2024's
(Z)-3-hexen-1-ol on PPI — neither convertible without a protein molar mass. **[P]**

### (j) WHAT CHANGED

**New:** `src/protein_binding.py`; `data/lit/binding_constants.yml`;
`scripts/generators/generate_matrix_binding_mode_comparison.py`;
`results/validation/matrix_binding_mode_comparison.{json,md}` (un-gitignored, with the reason
in place — and the filename must NEVER contain `refit`/`rederivation`, because
`holdout_guard` check 4 treats those as fit records and this file names every hold-out
benchmark); `tests/scientific/test_wave_s4_protein_binding.py` (18 tests).

**Edited:** `src/headspace.py` (one optional kwarg + the binding branch),
`src/benchmark_validation.py` (context resolution, the guard, two metadata fields),
`.gitignore`, `README.md`, `AUDIT.md`.

**NOT changed:** any observability factor, any marker yield, any barrier, any projection
constant, any benchmark, any hold-out bundle. `fitted_factors` is still the default and
`test_shipped_pea_prediction_is_unchanged_by_this_wave` pins the shipped 1125.278 ppb.

**Regeneration scope: none needed.** No shipped number moved, which is measured rather than
asserted — the `fitted_factors` column of the mode comparison reproduces
`external_validation_report.json`'s 93.6837x median and 3/8 coverage exactly, and 186 existing
matrix / headspace / benchmark-summary / honest-headline tests pass unchanged.

### (k) GATES + TESTS

| Gate | Result |
|---|---|
| `scripts/ci/holdout_guard.py` | **PASS** — 4/4; 16 bundles, 12 flagged, 34 fit records scanned |
| `scripts/ci/citation_gate.py` | **PASS** — **98** files (was 97), **1013** DOI-bearing fields (was 1002), **331** unique DOIs (was 323), 0 waivers. The 8 new DOIs are this wave's binding sources. |
| `scripts/ci/fit_target_gate.py` | **PASS** — 2/2. `binding_constants.yml` declares no fit target because there is none. |

Targeted tests only, per the owner directive (**no full-suite run**):
`tests/scientific/test_wave_s4_protein_binding.py` **18 passed**; the regression subset
`-k "matrix or headspace or honest_headline or benchmark_summary or observability or
external_validation or uncertainty"` **186 passed, 1 xfailed, 0 failed**. No existing test
was re-pinned, because no shipped number moved.

The new tests deliberately do NOT pin the mode-vs-mode scores — pinning a number the wave
just produced is the circularity Rounds 1-3 removed, and Wave U set that precedent. They pin
the inertness of the default path, the arithmetic and unit handling of the model, the
provenance discipline of the data file, and that the double-count guard fires.

### (l) [P] CARRIED FORWARD

- [P] **The absolute-scale deficit is in the marker yields, not in observability**, and now
  there is a unit argument for it: an observability factor cannot exceed 1 and two shipped
  ones are 4.32 and 9.54. `MATRIX_BENCHMARK_BASE_MARKER_YIELDS` (0.205 / 0.502 / 0.063 /
  0.150) is `260 / 638 / 80 / 150` over a common scale, i.e. still built from the retired
  values. This is the volatile-budget workstream's highest-value target.
- [P] **Adopt, keep or retire the fitted observability factors.** The measured comparison is
  on the record: they win in-panel by 5.9x in median fold (1.0004x vs 5.92x) and lose to *doing nothing*
  out of sample. This wave does not make that call; it makes it decidable.
- [P] **The 1-hexanol lane still has no anchor of any kind** (see (i)). Retiring it now has a
  measured consequence (1117x -> 492x on the only scored row) rather than a guess.
- [P] **Buy or ILL three papers.** `10.1016/j.foodchem.2022.133044` (Bi et al. 2022) is the
  only identified primary source with a hexanal x pea-protein binding constant and is closed
  on every open route tried (Unpaywall, OpenAlex, Semantic Scholar, Europe PMC, fatcat, IA
  Scholar); it is quoted second-hand by Wongprasert 2024 as **K ~ 684.46 M^-1, n = 4.58 at
  37 C**, with no conditions given. Also `10.1021/jf00006a003` (O'Keefe 1991 JAFC, to
  disambiguate the 160/320 kDa sentence) and `10.1021/jf00108a037` (Damodaran & Kinsella
  1981, quoted second-hand by two independent theses).
- [P] **A protein molar mass would unlock four more datasets.** Barallat-Perez 2025, Liu 2026,
  O'Keefe 1988 and Wongprasert 2024 all report K in M^-1. A REPORTED (not modelled)
  consistency calculation: matching Snel's per-gram `a_p` for soy 2-octanone against
  Barallat-Perez 2025's nK for the same compound implies an effective protein molar mass of
  ~53 kDa, which is physically sensible for a soy subunit. That is a unit reconciliation
  between two literature datasets, **not used anywhere in the model**, and it should not be
  adopted without a source that states the mass.
- [P] **No source anywhere gives 2-pentylfuran binding to any protein.** Both retrieval
  sweeps returned zero. Its `f_free` in this model is the least supported (Pow is XLogP3 and
  the class transfer is from ketones).
- [P] **Barallat-Perez 2023 and the a_p route disagree on the SIGN of the pea-vs-soy hexanal
  binding difference.** The paper says soy binds hexanal more than pea (qualitative, no
  number); `a_p` says pea binds more (25e-5 vs 16e-5); and the shipped registry says soy
  RELEASES 2.2098x more than pea. Two open-access sources contradict each other and neither
  prints a hexanal number. This is the sharpest remaining evidence gap in the lane.

---

## Wave S5 — the comparative front door, a generated model card, and nine fewer documents (2026-08-28)

**Mandate:** build one user-facing entry point around the product thesis the owner approved —
*comparisons cancel systematic error; the model's measured skill is ordinal* — generate the
validity domain instead of writing it, and finish with **fewer** documents than we started.
No new science, no re-tuning. All three constraints held; the file count is **−13**.

### (a) THE PREMISE, RESTATED FROM THE EVIDENCE RATHER THAN ASSERTED

| lane | out-of-sample absolute accuracy | ordinal accuracy |
|---|---|---|
| free precursor | median **6.04×**, worst 52.6×, 12/21 within 10× | 8/12 within-point orderings, **3/6** series directions |
| protein matrix | median **67.4–93.7×** across all three observability modes | — |
| genuine extrapolation | **1/5** inside the 90% CI | — |
| directional panel | — | **21/29** independent, **17/19** excluding pH/aw |

Nothing above is new. The wave's contribution is that a user now meets those numbers at the
point of use rather than in a README they may not read.

### (b) THE STALE-ARTIFACT FINDING — the reason the card is generated, not written

`docs/validation/directional_accuracy_report.md` was still publishing **20/29** as its
post-fix headline. The shipped tree scores **21/29**: Wave T4's Heyns-consistency fix flipped
SUG-12 (HMF, fructose vs glucose) MISS → OK, carrying `sugar_identity` 7/8 → 8/8, and Waves S3
and S4 both re-scored the panel byte-identically at 21/29. The ledger recorded it; the
artifact did not. **Three waves published a number one point below the truth, in the file
whose entire purpose is to state that number.**

A second, sharper version of the same defect was found while wiring the tags: the report's
**category table and its bucket table have different denominators.** The categories sum to
27/35 — the *screening* bucket (29 independent + 6 overlap) — so deriving "excluding pH and
aw" by summing categories yields 23/25 (92%), where the headline-comparable figure is
**17/19 (89%)**. Both are true; mixing them inflates the score by three points, and the first
draft of the model card did exactly that. Fixed by adding two explicitly-labelled bucket rows
to the report and having the generator **read** them rather than derive them, with the
arithmetic written out in place. All ten pH/aw rows are independent, which is why 19 + 10 = 29
exactly and why the six overlap rows all sit outside pH and aw.

`docs/validation/directional_accuracy_report.md` gains a `## CURRENT STANDING` section — the
one place stating the live counts, machine-read by the CLI and the card — plus a pointer at
§2 so no reader stops at the Wave S2 numbers. **Everything above that section is left
verbatim; the pre-fix numbers are the evidence for the delta.**

### (c) THE COMPARATIVE INTERFACE — `python scripts/maillard.py`

Three verbs, thin orchestration over existing machinery:

| verb | leads with | reliability |
|---|---|---|
| `compare` | per-compound **A/B ratios**, dominant (rate-limiting) pathway per side | per-axis tag from the panel, weakest axis governs |
| `predict` | a **range**, not a point; model-card caveats inline | lane tag per compound |
| `rank-experiments` | the VoI queue, promoted to a first-class verb | — |

Absolute ppb appear only behind `--absolute`, and the caveat ships **in the same block** as
the numbers, never as a footnote. `--json` emits the payload; `--top N` trims.

**The reliability tags are read from the artifact at runtime, not compiled in.** Re-score the
panel, edit the CURRENT STANDING table, and the tool's advice moves — there is no second copy
of the numbers. A missing or unparseable report **raises**; it does not fall back to a baked-in
default, because a silent default is precisely how a number outlives its evidence.

Thresholds, stated in the card and applied without judgement: `trust` ≥ 80% on ≥ 3 claims,
`caution` ≥ 60% or too few claims, `do-not-use` < 60% **or unmeasured**. Under those, a sugar
swap reports `trust (8/8)`, a pH sweep `do-not-use (4/7)` and a moisture change
`do-not-use (0/3)` — the conclusion §A6 reaches in prose, now enforced at the point of use.
Bundling a pH change with a sugar change does **not** launder the pH miss: the weakest axis
governs, and that is a test.

**What the ratio interval is, and is not.** `A_p5/B_p95 .. A_p95/B_p5` — the ordinary
independent-error combination of two existing envelopes, labelled in output as the
conservative bound it is. The comparative thesis says the two arms' errors are strongly
correlated and the true interval is far tighter. **Nobody has measured that correlation, so
this wave models none.** Stating the bound and naming it as too wide is the honest position;
inventing a correlation coefficient to tighten it would have been new science.

**One plumbing change to `src/`, no science:** `MaillardPipeline.last_route_traces`, a
side-channel exactly like the existing `last_skipped_formulations`, exposing the
`debug_paths` / `debug_channel_flux` payload `predict_from_steps` already returned and
`evaluate_all` already discarded. Not a `FormulationResult` field, so no debug payload enters
any serialized report. Measured inert: the 281-test regression subset is green and the
determinism test still passes.

**A presentation defect found and fixed:** `predicted_ppb` carries **three** aliases per
compound (canonical SMILES, species name, display label) for downstream benchmark matching.
Read naively, the table printed every compound three times. Rows are now built from
`result.targets`, which is already one per compound. The same aliasing is why the sulfur and
lipid caveats match on substrings — a caveat that fires on only one spelling is a caveat that
silently does not fire.

**Also surfaced rather than smoothed over:** `predict` prints the run-level tier, which for a
cys/ribose system reads `high`. That is the *scope* vocabulary, not a validation claim, and
the README says there is no high trust tier. The renderer now says so on the line below,
naming 0/14 strict-ready — one honest sentence beats renaming a field the rest of the
codebase uses.

### (d) THE MODEL CARD — `scripts/generators/generate_model_card.py` → README.md

Reads the directional report, both hold-out artifacts, the mode comparison, the benchmark
panel (**recomputed live**, because `benchmark_summary.json` is gitignored and a card that
read it would print "artifact absent" on every fresh clone — the self-excusing-skip pattern
Wave J2 removed), the `no_verifiable_source` census (recounted), the sulfur benchmark's anchor
status (**checked, not asserted**) and the three gates (**actually run**). Emits a claim-type ×
system-class table with trust/caution/do-not-use verdicts, spliced between
`<!-- BEGIN GENERATED: model-card -->` markers. **No new standalone document**; the
machine-readable twin lands in gitignored `results/validation/model_card.json`.

The repo had **no** marker precedent — `README.md` contained zero HTML comments and no
generator has ever written it (`generate_report_visual_examples.py` only *prints a reminder*
to copy by hand). This wave establishes the convention. `--check` makes drift a failing
command; a test asserts the committed card equals what the artifacts currently say.

**The three honest headline sentences are computed, and one of them is verified rather than
recited:** the sulfur claim only renders if `cys_ribose_140C_Hofmann1998.json` actually shows
every `measured_volatiles` entry carrying `value_status: no_verifiable_source`. If a real
anchor ever lands, the card stops making the claim by itself.

**A census contradiction was caught before it shipped.** The recount finds **166** records
carrying `no_verifiable_source` in a status field; README's test-pinned figure is **120**.
Both are right under different definitions — 120 is the `source_status` key specifically, which
the recount reproduces **exactly**. The card now leads with 120, states the 46-record
remainder and its keys, and points the 98-numeric / 80-runtime split at the guard test that
owns it. Publishing 166 beside a pinned 120 would have read as drift in a document whose
purpose is to prevent drift.

**A negative result worth recording:** no artifact state can make the card say `trust` about an
absolute concentration. That is a test, not a convention.

### (e) CONSOLIDATION — 30 → 17 documents, 9157 → 7387 lines

| document | disposition | lines |
|---|---|---:|
| `docs/research/archives/` ×5 digests | **DELETED** — md5-identical to `data/Gemini_Deep_Research/` | −1529 |
| `docs/research/archives/README.md` | folded into `data/Gemini_Deep_Research/README.md`, **DELETED** | −66 |
| `docs/Plant Protein Thermal Processing Research.md` | **DELETED** — a *sixth* md5-identical duplicate, zero references, and the only one with **no warning header at all** | −173 |
| `docs/architecture.md` | folded into README (corrected), **DELETED** | −195 |
| `docs/reference/SMIRKS_SYSTEM.md` | folded into README (corrected), **DELETED** | −37 |
| `docs/reference/COMMAND_REFERENCE.md` | folded into QUICKSTART appendix, **DELETED** | −85 |
| `docs/guides/CURATING_LIPID_OXIDATION_ANCHOR.md` | folded into VALIDATION_CONTRACT §3E-bis, **DELETED** | −79 |
| `docs/notebooks/2_Food_Scientist_Walkthrough.md` | **DELETED** — a changelog for a notebook nothing references | −87 |
| `Max_ref.azure.txt` | **DELETED** — a tracked 0-byte file at the repo root | −0 |
| README.md | +ARCHITECTURE, +model card, +CLI Step 0 | +204 |
| QUICKSTART.md | +command reference, +CLI section | +125 |
| VALIDATION_CONTRACT.md | +anchor-curation methodology | +76 |
| directional_accuracy_report.md | +CURRENT STANDING | +66 |
| **NET** | **−13 documents** | **−1770** |

**KEPT, deliberately, against the count.** `docs/guides/GLOSSARY.md` (71 lines) was a
fold candidate and is not folded. It is the corrected post-audit version, it is the only place
the three colliding tier vocabularies are disentangled at length, and pouring 71 dense lines
into a 1300-line README makes the README worse. *The directive is fewer documents, not fewer
good documents;* there were enough genuine husks to make the number without deleting a live
one. Recorded here so the decision is visible rather than silently taken.

**The duplicate deletion is a provenance question, and it was treated as one.** The archive
README argued the five copies were "the evidence trail" for findings that trace a shipped
constant to a digest and nothing else. That argument protects the **content**, and the content
is untouched: every file was byte-identical (`cmp`-verified, not assumed) to a copy under
`data/Gemini_Deep_Research/`, which carries the stronger warning README. Deleting a second copy
destroys no evidence; keeping an unwarned copy under an authoritative-looking `docs/` path was
itself the defect Wave T3 opened, and that file's own "Preferred long-term disposition" section
prescribed this consolidation. Its unique finding — the dangling `pathways.md` — was folded
into the surviving README **before** the delete.

**Eight citation strings repaired, two of which were already broken** before this wave:
`data/lit/computational_priors.json` cited `docs/Maillard_Plant_based.md` and
`data/lit/matrix_family_coverage_registry.json` cited `docs/Elicit - … Report.md`; **neither
path has ever existed** (both files lived under `docs/research/archives/`). Repointing them is a
fix, not a move. **`pathways.md` was NOT repointed** — it exists nowhere, and the three headers
naming it supply the odour thresholds that are the denominator of every OAV in
`src/sensory.py`. Its name is kept and flagged in place in all three files. Substituting a
plausible source for an absent one is the defect this remediation exists to remove. Carried
`[P]`, unchanged.

`citation_gate.DIGEST_SOURCE_RE` keeps its `docs/research/archives` alternative even though the
directory is gone, with the reason written in place: a stale citation string still naming that
path must keep tripping the detector rather than reading as clean.

**Pre-audit claims corrected in the folds, not carried across.** `docs/architecture.md` opened
its trust surface with "**High trust — use freely**" and advertised the retired
`bounded_calibration` / `transferred_literature` / `surrogate_family` vocabulary;
`SMIRKS_SYSTEM.md` called two templates "Validated vs. benchmarks" at 0/14 strict-ready. Both
claims were dropped, not relocated. What survived: the **two-lane design** (`free_precursor` vs
`matrix_only`, and why three network waves left all eight matrix hold-out points bit-identical),
`trunk_kinetics` and why it deliberately ships unwired, the benchmarks→gates→artifacts→README
data flow, and the SMIRKS Tier A/B provenance with the MW-300 Da cap stated as the **scientific
limitation** it is (it is why 13-HPODE at 312 Da is pruned and the lipid chain cannot reach
hexanal).

**Two stale links fixed in passing**, both pre-existing: `SCIENTIFIC_REFERENCE.md` pointed four
links at `docs/use_cases/`, a directory that has never existed (the files are in
`docs/protocols/`), and `data/lit/README.md` carried an absolute
`file:///Users/pabloantoniomorenocasares/...` URI. `CONTRIBUTING.md`'s `docs/` map omitted four
real subtrees and named three files that no longer exist; it also asserted "No Python execution
code" under `data/` while `data/reactions/curated_pathways.py` sits there. Corrected in place.

### (f) GATES + TESTS

| Gate | Result |
|---|---|
| `scripts/ci/holdout_guard.py` | **PASS** — 4/4; 16 bundles, 12 flagged, 34 fit records scanned |
| `scripts/ci/citation_gate.py` | **PASS** — **92** files (was 98), **909** DOI-bearing fields (was 1013), **287** unique DOIs (was 331), 0 waivers. The drop is entirely the six deleted duplicate digests; no DOI was removed from any *data* file, and the surviving copies carry the same DOIs under `data/Gemini_Deep_Research/`. |
| `scripts/ci/fit_target_gate.py` | **PASS** — 2/2 |

Targeted only, per the owner directive (**no full-suite run**):

* `tests/unit/test_comparative_cli_2026_08.py` — **39 passed** (new: CLI shape per verb,
  rendering honesty, reliability wiring, model-card generation and staleness)
* affected docs/front-door tests — `test_honest_headline_guards.py`,
  `test_front_door_ux_2026_08.py`, `test_wave_i_chemistry.py`: **121 passed**
* regression subset `-k "matrix or headspace or benchmark or observability or
  external_validation or uncertainty or honest_headline"`: **281 passed, 1 xfailed, 0 failed**

The new tests deliberately pin **no predicted ppb, ratio or fold error** — pinning a number
this wave produced is the circularity Rounds 1–3 removed. They pin behaviour: that the tags
follow the artifact (verified by pointing the loader at a synthetic report and watching the
verdicts invert), that a missing report raises rather than defaults, that absolutes never reach
a user without their caveat, that unmeasured reads as `do-not-use`, and that no artifact state
can produce a `trust` verdict on an absolute concentration.

**One packaging defect caught at the end.** `data/cli_examples/compare_ribose_vs_glucose.yml`
— the worked example the CLI smoke test runs against — was caught by `.gitignore:33`'s
blanket `data/*`. Untracked, its test would have passed locally and failed on a fresh clone:
the self-excusing-skip pattern Wave J2 removed. Un-ignored with the reason written in place.

### (g) [P] CARRIED FORWARD

- [P] **`pathways.md` still does not exist** and is still cited by three runtime-critical data
  files. Unchanged by this wave by design; re-sourcing per-compound odour thresholds is a
  retrieval job, not a documentation one.
- [P] **`docs/notebooks/results/maillard_results.db`** is a generated SQLite artifact tracked
  inside `docs/`, contrary to CONTRIBUTING's own "`results/` is gitignored" rule. Left in place:
  a notebook may read it, and deleting tracked binary data on a documentation pass is the wrong
  wave for that call.
- [P] **Neither `docs/notebooks/*.ipynb` is referenced by anything**, and their committed
  outputs may embed pre-audit numbers. No gate scans notebook outputs.
- [P] **The ratio-interval correlation is unmeasured.** The comparative thesis is that errors
  cancel between arms; the CLI can only publish the bound that assumes they do not. A
  measurement of that correlation — e.g. paired hold-out arms scored as ratios — would be the
  single highest-value addition to the interface, and would convert a caveat into a number.
- [P] **`docs/reference/SCIENTIFIC_REFERENCE.md` is date-stale** ("as of 2026-03-18"): its
  "well supported" and numeric-anchor tables predate the entire citation audit and have never
  been reconciled against it. It is referenced by path from `src/reporting.py:1102`, so it
  cannot simply be deleted; it needs a content pass.
- [P] **`docs/slr_benchmark_evaluation.md` carries a dead branch reference** (`real_tests`) and
  a March date header, while two `data/lit/` registries name it as their `source_document`.

---

## Wave W — the ILL payoff: the sulfur branch gets real anchors, and fails them

**2026-08-28.** Five previously-paywalled papers were delivered by the owner into
`data/articles/` (gitignored; **no PDF is committed**). This wave read all five in full and
turned them into data. It fitted **nothing** — no barrier, no observability factor, no
threshold, no prior was touched anywhere in the tree.

### The headline, stated the unflattering way round

**"THE SULFUR BRANCH HAS ZERO ABSOLUTE LITERATURE ANCHORS" IS NO LONGER TRUE. It now has
three, and the model fails all three.** That sentence had been the audit's signature finding
since Wave S2c; it is corrected in README.md, AUDIT.md, VALIDATION_CONTRACT.md,
`src/barrier_constants.py`, `src/comparative_cli.py`, `src/model_card.py`,
`src/curated_pathways.py`, `data/benchmarks/maillard_validation_benchmarks.md` and the two
tests that pinned it. The three anchors are the pH-5.0 aqueous rows of Hofmann & Schieberle
1998 Table 1, and the model misses them by **12.27×, 29.58× and 14.46×** (six comparisons,
mean 0.92 dex). *The branch went from unanchored to anchored-and-measurably-wrong, which is
the direction of travel this audit wanted.*

### (a) The anchor, and what the primary source settled

`hofmann1998_ribose_cysteine_145C_20min_pH5` — **FFT 121 ppb, MFT 198 ppb**, from cysteine
3.3 mmol + ribose 10.0 mmol in **100 mL** of 0.5 mol/L phosphate at pH 5.0, **145 °C / 20 min**
in a laboratory autoclave, quantified by SIDA against [²H₂]-FFT and [²H₃]-MFT. Native ppb is
legitimate here and this is the whole difference from the retired file: **the volume is
printed in the table footnote**, so µg-per-100-mL × 10 = µg/L is a single stated
multiplication with no assumed concentration. The retired benchmark's ppb needed a
mol%→ppb conversion resting on an unattested 10 mM basis — a free multiplicative parameter.

**Every table was independently re-read from the PDF by this wave, not taken on trust from
the orchestrator's transcription.** All ten matched.

Three things the primary source settled that no amount of secondary retrieval could:

1. **342 and 200 ppb appear nowhere in the paper.** Wave S2b's forensic conclusion — that
   they were a repo-internal derivation from `maillard_validation_benchmarks.md` §1.3 — is
   **confirmed from the source**. The two files are cross-linked in both directions;
   `cys_ribose_140C_Hofmann1998.json` is kept, unedited, as the provenance record, and its
   numbers are not resurrected, copied or reconciled.
2. **§1.3's invented band was not merely unsourced, it was wrong by an order of magnitude,
   in the optimistic direction.** Its `~0.02–0.05 mol %` MFT for ribose+cysteine, checked
   against the real Table 1 (19.8 µg from 10.0 mmol): **0.00173 mol %**, about **12× below
   the bottom of the band.** Recorded in `maillard_validation_benchmarks.md` itself.
3. **The conditions block was not even a faithful copy of the protocol it was borrowed
   from.** Mottram & Nobrega's real system is 140 °C / 30 min at **50 mM** each; the retired
   file declared 10 mM. Wave S2b could not check that; now it is checked, and it is wrong too.

**The contract, derived from the measurement and from nothing else.** Table 1 footnote b:
*"Data are mean values of triplicates and differed by not more than 10%."* Every other table
in the paper references it. 10% ⇒ `max_ratio` **1.10**, `mean_abs_log10_error` **0.0414**
(= log₁₀ 1.10). It is set **at** the measurement's precision, not above it: the failure mode
being cleaned up is the opposite one — the retired 1.45× / 0.09 dex was ~1.7× *tighter* than
its own invented source's spread. Note the shipped contract is **stricter** than the global
free-precursor default (1.5 / 0.10), which was available and would have been looser. Nothing
was widened so a row would pass. The model fails both.

### (b) The fit/hold-out split — decided and written into the files BEFORE ingestion

**The split is on an axis, not on convenience: the scored panel is the pH-5.0 slice; the
hold-out is everything off pH 5.0, plus one sugar held back.**

| | rows | why |
|---|---|---|
| **PANEL** `data/benchmarks/` | ribose, glucose, fructose @ pH 5.0 | the classic condition, and the one the rest of the panel already sits at |
| **HOLD-OUT** `external_validation/maillard_path/` | xylose @ pH 5.0; ribose & glucose @ pH 3.0 and pH 7.0 | xylose is the one executable sugar the panel never sees; the pH axis is the pre-registration target of the pending pH/aw work |

The load-bearing half is the pH series. Table 2 contains a **sign-crossing** claim — from
ribose both odorants rise monotonically as pH falls 7→5→3, while from glucose and rhamnose
both go through a **maximum at pH 5.0** — and that is precisely the acid test of any future
pH term. Putting it in the panel would make it fittable; putting it in the hold-out
pre-registers it. **No row is on both sides. Nothing was fitted to anything.**
`holdout_guard` passes (21 bundles, 17 maillard-path).

**The dry-heat rows are deliberately NOT ingested as absolute rows anywhere.** The paper
states no water activity and no system mass for the silica-gel system, so both the engine's
`water_activity` input and the µg→ppb denominator would have to be **invented** — the same
class of move that produced 342/200. They are ingested as *ordinal* claims instead, which
need no such number.

### (c) Model vs measurement — every new row, scored once and frozen

| | benchmark | analyte | measured | predicted | fold | \|dex\| |
|---|---|---|---:|---:|---:|---:|
| PANEL | ribose pH 5.0 | FFT | 121 | 1485.08 | **12.27×** | 1.089 |
| PANEL | ribose pH 5.0 | MFT | 198 | 468.58 | 2.37× | 0.374 |
| PANEL | glucose pH 5.0 | FFT | 28 | 828.35 | **29.58×** | 1.471 |
| PANEL | glucose pH 5.0 | MFT | 19 | 90.23 | 4.75× | 0.677 |
| PANEL | fructose pH 5.0 | FFT | 32 | 5.39 | 0.17× | 0.774 |
| PANEL | fructose pH 5.0 | MFT | 25 | 1.73 | **0.07×** | 1.160 |
| HOLD-OUT | xylose pH 5.0 | FFT | 96 | 1485.08 | 15.47× | 1.189 |
| HOLD-OUT | xylose pH 5.0 | MFT | 143 | 468.58 | 3.28× | 0.515 |
| HOLD-OUT | ribose pH 3.0 | FFT | 229 | 1662.91 | 7.26× | 0.861 |
| HOLD-OUT | ribose pH 3.0 | MFT | 553 | 2145.61 | 3.88× | 0.589 |
| HOLD-OUT | ribose pH 7.0 | FFT | 12 | 641.02 | 53.42× | 1.728 |
| HOLD-OUT | ribose pH 7.0 | MFT | 25 | 161.32 | 6.45× | 0.810 |
| HOLD-OUT | glucose pH 3.0 | FFT | 7 | 951.23 | 135.89× | 2.133 |
| HOLD-OUT | glucose pH 3.0 | MFT | 3 | 1519.27 | **506.42×** | 2.705 |
| HOLD-OUT | glucose pH 7.0 | FFT | 6 | 343.39 | 57.23× | 1.758 |
| HOLD-OUT | glucose pH 7.0 | MFT | 4 | 22.85 | 5.71× | 0.757 |

**Mean |log₁₀| — panel 0.924 dex (6 rows), hold-out 1.304 dex (10 rows).** The hold-out is
0.38 dex worse than the panel on rows from the *same paper, same table, same lab, same day*.
That gap is not out-of-distribution generalisation failure; it is the model degrading as it
moves off pH 5.0.

#### The pentose ≫ hexose test — the first same-conditions measurement this claim has ever been scored against

| ratio | measured (Table 1) | model | verdict |
|---|---:|---:|---|
| **MFT ribose/glucose** | **10.42×** | **5.19×** | model has **half** the real pentose advantage |
| **FFT ribose/glucose** | **4.32×** | **1.79×** | model has **41%** of it |
| MFT ribose/fructose | 7.92× | **270.99×** | 34× too large |
| FFT ribose/fructose | 3.78× | **275.64×** | 73× too large |
| MFT ribose/xylose | 1.38× | **1.0000×** | model cannot distinguish two pentoses at all |
| FFT ribose/xylose | 1.26× | **1.0000×** | ditto |

The repo's own pentose≫hexose figure has ranged **3.4×–18.3×** across waves; the measurement
is **10.42×** (MFT) and the shipped model says **5.19×**. So the claim is real, the model has
the sign right, and it **under-states the effect by about 2×**. That is the most useful thing
in this wave: a number that was previously only bounded by which wave you asked.

**And two failures that were invisible before:**

- **The model invents a 154× gap between two hexoses that measure within 1.3× of each other.**
  Fructose vs glucose: measured FFT 3.2 vs 2.8 µg, MFT 2.5 vs 1.9 µg. Model: glucose ahead by
  52× (MFT) and 154× (FFT), with the **sign wrong** as well. Being good at large real gaps and
  inventing large fake ones is the same defect seen from two sides.
- **Ribose and xylose are byte-identical to six significant figures** — the network reaches
  both through one pentose lane. The hold-out row is structurally guaranteed to fail on the
  ordinal axis, and that is the finding, not a defect in the row.

#### The sign-crossing pH claim (hold-out): the model gets the pentose right and the hexose wrong

| | measured pH 3 / 5 / 7 | shape | model | shape | |
|---|---|---|---|---|---|
| ribose FFT | 229 / 121 / 12 | monotonic ↓pH ⇒ ↑ | 1663 / 1485 / 641 | monotonic | ✅ |
| ribose MFT | 553 / 198 / 25 | monotonic | 2146 / 469 / 161 | monotonic | ✅ |
| glucose FFT | 7 / **28** / 6 | **max at pH 5** | 951 / 828 / 343 | monotonic | ❌ |
| glucose MFT | 3 / **19** / 4 | **max at pH 5** | 1519 / 90 / 23 | monotonic | ❌ |

The model applies one pH shape to both sugars, which is exactly what a family-level pH term
must do — and exactly what the measurement forbids. **This is now pre-registered, out of
sample, and cannot be fitted to.**

### (d) The step-level gap list — what the engine cannot execute, and what it would need

Tables 3, 4, 6, 7, 8 and 10 are the repo's **first constraints on individual steps rather
than end-to-end lumps**, and they are the reason this paper was worth the ILL. **None of them
is executable today.** `src/precursor_resolver.py` resolves **19** species
(`data/species/precursors.yml`); every step-level precursor in the paper is missing:

| precursor | Hofmann tables | status | what the row would constrain |
|---|---|---|---|
| hydroxyacetaldehyde (glycolaldehyde) | 6, 8, 10 | **MISSING** | the C2 half of Wave P's C2+C3 lane |
| 2-oxopropanal (methylglyoxal/pyruvaldehyde) | 6, 7 | **MISSING** | the C3 half |
| mercapto-2-propanone | 7, 8, 9, 10 | **MISSING** | Wave P's lane, measured end to end |
| norfuraneol (NF) | 4, 5, 10 | **MISSING** | the route Wave N retired — Table 5 is the evidence |
| 3-deoxyribosulose | 3 | **MISSING** | the deoxyosone branch point |
| furan-2-aldehyde (furfural) as a *precursor* | 3 | **MISSING** | FFT's dominant route (60× ribose) |
| hydrogen sulfide (H₂S) | 3, 4, 7 | **MISSING** | the sulfur donor in every step-level row |
| rhamnose, glucose-6-phosphate, maltose | 1, 2, 5 | **MISSING** | 3 more Table 1 rows |
| ribose-5-phosphate, IMP | *(Mottram)* | **MISSING** | the meat-relevant precursor pool |
| thiamine | 10 | present | but Table 10's row is single-precursor |

**What would be needed, precisely:** adding a species to `precursors.yml` is *not* sufficient
and would be actively harmful — the resolver would then return a confidently wrong answer for
a precursor the reaction network has no route from. Each row needs the species **and** the
step: e.g. ribose-5-phosphate needs a dephosphorylation step (the paper's own mechanism), and
furan-2-aldehyde-as-precursor needs the H₂S addition step that Table 3 measures.

**Two things the paper hands us for free, recorded now so a later wave does not have to
re-derive them:**

- **The mol% basis is documented.** mol% = mol product / mol limiting precursor × 100,
  verified against Table 4 (ribose 15.1 µg MFT / 114.17 g·mol⁻¹ = 1.323e-7 mol; per 1 mmol =
  0.0132% ≈ the printed 0.01%). Any future mol% target now converts honestly.
- **Table 5 is independent support for Wave N's route correction, from the very paper whose
  misattributed values funded the norfuraneol barrier fit.** Ribose makes **54 530 µg NF** and
  only **19.8 µg MFT** — a factor of ~2750 — while glucose's ratio is ~7. The authors' own
  conclusion: *"MFT formation might not run exclusively via NF as the key intermediate."*
- Table 10 ranks the routes: hydroxyacetaldehyde/mercapto-2-propanone **268.1 µg (0.24 mol%)**
  > NF/H₂S 211.2 > NF/cysteine 50.8 > **thiamin 8.2**. **Wave P was right to add the C2+C3
  lane** — it is the single most effective measured MFT route, and thiamin is 33× weaker.

### (e) What the other four papers unlocked — every old → new

#### Bi et al. 2020 (`10.1021/acs.jafc.9b07711`) — the two-branch question, answered

Wave I flagged that the hold-out hexanal value 1260 ppb was **exactly** OAV 280 × the repo's
own unsourced 4.5 ppb threshold, and could not tell from inside the repo whether the paper
reported the concentration or the OAV. **The paper answers it: the concentration is the
measured quantity and the OAV is derived from it.** Table 4's footnote h defines OAV as the
*"ratio of concentration to odor threshold"*; the concentration column carries a standard
deviation and the OAV column does not. The repo was running the paper's own arithmetic
backwards and landing on the right number by construction.

| | old | new | fold change |
|---|---|---|---|
| raw pea hexanal | 1260.0 ppb, `derived_from_oav_and_repo_threshold` | **1260.0 ppb**, `primary_source_table_value_verified`, ±114 (9.05%) | **1.000×** |
| roasted pea hexanal | 324.0 ppb, same | **324.0 ppb**, same, ±34.3 (10.59%) | **1.000×** |

**The number does not move. This is a provenance correction and an uncertainty addition, and
it is reported as that rather than dressed up as a fix.** What changes is what the value is
entitled to be called: it was "a consistency check against our own threshold"; it is now an
external measurement with a stated replicate spread. **The coincidence is explained**: Bi et
al.'s own threshold column prints **4.5 µg/kg** for hexanal, attributed to van Gemert — the
same compilation the repo uses, which is why the division came out exact. That corroborates
the repo's 4.5 ppb ODT's provenance; it is **not** an independent measurement of it.

Quantification class recorded honestly: external calibration curve against
2,4,6-trimethylpyridine, R² 0.997, triplicate — **not** isotope dilution, so weaker than the
Hofmann and Yiltirak sulfur anchors, and labelled as such rather than levelled up.

#### Mottram & Nobrega 2002 (`10.1021/jf0200826`) — ORDINAL, permanently

**Wave S2b predicted "Tenax-relative at ~85% confidence". CONFIRMED, and the paper says it
itself:** *"This allowed comparison of the relative contributions the volatiles made to the
headspaces of the different systems **but did not provide absolute concentrations in the
aqueous solutions**."* Table 1's footnote: *"Approximate quantities in headspace (ng/mmol of
sugar)"*. **It can never become an absolute benchmark here, at any future date** — the
quantity it reports is not the quantity the model predicts. Six ordinal claims added to
`docs/validation/directional_claims_panel.yml`, where the panel is weakest.

#### Brands & van Boekel 2001 (`10.1021/jf001430b`) — 39 trajectory series

`data/lit/timeseries/brands_sugar_casein_120C_pH68.yml` held **endpoints only** and said so:
*"Reading points off those figures would be eyeball estimation at low resolution and was NOT
attempted. Those numbers are NOT in this repository."* They are now — **not** from the thesis
scan, but from the journal paper's own figures. **The paper contains ZERO tables**, so
`data_quality: table` is unachievable from this source by anyone, ever; every series is
labelled `figure_readoff` under a new `source_2001` block kept separate from the thesis data.

**Four independent verification checks, all passed and all recorded:** (1) the digitized
"acids by HPLC" series from Fig. 3 equals formic + acetic digitized separately from Figs. 1/2
to within 2% at every point; (2) the digitized pH endpoints reproduce the paper's prose pH
drops (0.3 and 0.4 units) to 0.03; (3) Fig. 4's two y-axes stand in exactly the 1:2 ratio the
paper's own ε = 500 implies; (4) **this wave re-read Figure 7 and Figure 1's left panel by
eye from upsampled renders, without reference to the digitizer** — Fig. 7 agreed on all 18
values to one least-significant digit, Fig. 1's glucose series on all 8 to ~1 mmol/L.

Four superimposed points are carried `ambiguous: true` rather than resolved; one point (Fig. 6
fructose formic acid at t=15) is **absent from the paper** and is omitted, not interpolated.
Stated plainly: the ±1% digitization error **is not the experiment's error bar** — the paper
reports no error bars at all.

Three things the file now holds that it did not: the **Amadori intermediate maximum**
(rises to ~1.1 mmol/L at t=10–15, then falls — a direct observation of an intermediate being
consumed, which no endpoint can show); the **aldose/ketose contrast** (Amadori flat at zero
for the whole fructose run against a peak for glucose); and the **no-protein caramelisation
control** (browning is ~6× protein-dependent for glucose but only ~2× for fructose).
`reports_no_kinetic_parameters` is recorded loudly: the paper fits **no** rate constants and
says so four times — if any constant here is ever attributed to this DOI, it is misattributed.

#### Bi et al. 2022 (`10.1016/j.foodchem.2022.133044`) — identity verified, and the answer is "no"

Identity confirmed from the PDF: *"Non-covalent interactions of selected flavors with pea
protein"*, Food Chemistry **389** (2022) 133044. **The hexanal × pea protein constant this
campaign most wanted is `K = 684.46 ± 109.43 M⁻¹`** (Table 1, 37 °C, pH 7.6) — **and it is
NOT usable**, for exactly the reason `binding_constants.yml` already documented for
`wongprasert2024_ppi_klotz`: K is defined per **mole** of protein and **the paper never states
a protein molar mass**. Promoting it would put a free multiplicative parameter under the most
load-bearing constant in the matrix lane. Recorded in
`not_usable_without_protein_molar_mass`, with the numbers, so a later wave can promote it.

What **is** usable is the PRV retention percentage, measured at a stated g/L loading and
needing no molar mass: **hexanal 71.72% bound at 10 g/L pea protein**, pH 7.6, 37 °C
(`bi2022_hexanal_pea_protein_prv`), plus (Z)-2-penten-1-ol 29.29% and (E)-2-octenal 79.31%.
It sits at the real in-pea concentration (~12.6 µM), far inside the dilute-ligand limit the
Langmuir form requires.

**Did the binding mode's hold-out score improve? NO — and it could not have.**

| mode | median fold error | before Wave W |
|---|---:|---:|
| fitted_factors | 93.68× | 93.68× |
| unit_observability | 67.42× | 67.42× |
| binding_physics | 68.18× | 68.18× |

**Unchanged to every digit.** Two reasons, both structural: the new record enters the
cross-check, not the prediction path (the shipped mode predicts from `a_p × Pow` alone, and
patching in a per-compound measured constant for hexanal *only* would make the mode a hybrid
of measured physics and per-compound patching — the exact pattern this audit removes); and
the Bi 2020 hold-out lanes are **raw pea flour**, for which no protein loading can be
asserted, so the binding mode applies no correction there at all.

**What the real constant did instead is falsify the binding model on the compound that
matters most.** The zero-parameter cross-check gained **its first hexanal row**:

> `bi2022_hexanal_pea_protein_prv` | hexanal | 10.00 g/L | measured **71.72%** | predicted
> **13.09%** | residual **−58.63 points**

That is **by far the largest residual in the table** (previous worst: +22.22). And the paper
explains it: *"hydrogen bonding was dominant in the non-covalent interactions between hexanal
and pea protein"*, in explicit contrast to the two other compounds it studied. **The shipped
binding physics is a purely hydrophobic partition model driven by Pow — so the one compound
the matrix lane leans on hardest is the one its functional form is least entitled to
describe.** Recorded as a mechanism-level caveat in `qualitative_only`; nothing was changed
for it, because no better form is specified well enough to implement. Two further caveats
recorded: the paper's own Ksv (14.205 L/mol) and Klotz K (684.46 M⁻¹) disagree **~48×** at the
same temperature in the same buffer and the authors never mention it; and binding percentage
**collapses** above 0.25 mM (a rise-then-fall a single-site Langmuir cannot produce at all).

### (f) Panel and directional aggregates — what moved, and which direction

| aggregate | before | after | why |
|---|---|---|---|
| Panel size | 14 | **17** | 3 Hofmann anchors |
| Strict-ready | 0/14 | **0/17** | all three new anchors FAIL |
| Predictive benchmarks w/o blocking gaps | 0/6 | **0/9** | denominator +50%, numerator still zero |
| MC panel | 11 benchmarks, 35 rows | **14, 41** | |
| Honest external-literature CI coverage | 0/3 @ **0.75 dex** | **4/9 @ 2.63 dex** | **see below — this is not an improvement** |
| Directional, strictly independent | 21/29 (72%) | **24/35 (69%)** | new rows scored 3/6 |
| …excluding pH and aw | 17/19 (89%) | **18/23 (78%)** | |
| `sugar_identity` | **8/8 "Trustworthy"** | **9/11** | first two misses ever |
| `ph` | 4/7 `do-not-use` | **6/9 `caution`** | both new pH rows AGREE |
| `ranking` | — | **0/1** | new category, one row, a miss |

**The CI-coverage move must not be quoted alone.** The old guard's message said: *"If this
rises, verify it is the model getting closer and not the interval getting wider — check
`median_ci_width_log10` in the same breath."* **It is the interval getting wider**: 0.746 →
2.627 dex, i.e. from a factor of 5.6 end-to-end to a factor of **~424**. Four of the six new
rows "cover" their measurement with intervals **2.90, 2.95, 3.04 and 3.74 dex** wide — an
interval spanning three orders of magnitude contains almost anything, so those hits are close
to vacuous. **No prior, threshold or interval was widened; the population changed, not the
model.** The two rows that *miss* (fructose FFT and MFT) are the informative half. The guard
now carries a two-sided width assertion so the vacuity cannot be quietly quoted away.

**The `ph` promotion is real but small, and came from evidence, not from a threshold.**
`CAUTION_MIN_RATE` is still 0.60 and `TRUST_MIN_RATE` still 0.80 — the test asserts both.
6/9 = 0.667 clears the floor by two rows, both from one source. `caution` here means "no
longer demonstrably worse than a coin", **not** "trustworthy"; the §8/§A6 licence still says
the model licenses no pH recommendation, and the re-pinned guard now asserts the invariant
that actually matters — **pH and water activity are never `trust`** — rather than one
particular tag.

**The `sugar_identity` demotion is the one to read.** It was 8/8 and labelled *Trustworthy*
without qualification. Every one of those eight rows was pentose-vs-hexose or thiol-vs-thiol
— comparisons across a large real gap. Given a hexose-vs-hexose comparison for the first
time, the model fails it by two orders of magnitude with the sign wrong.

### (g) Files

**New (8):** `data/benchmarks/hofmann1998_{ribose,glucose,fructose}_cysteine_145C_20min_pH5.json`;
`data/benchmarks/external_validation/maillard_path/mp_holdout_hofmann1998_{xylose_…_pH5,ribose_…_pH3,ribose_…_pH7,glucose_…_pH3,glucose_…_pH7}.json`.

**Corrected:** the two `external_validation_bi_2020_*_hexanal.json` bundles (provenance +
uncertainty, values unchanged, prior notes retained verbatim); `binding_constants.yml`
(+1 source, +3 percent-bound records, +1 not-usable Klotz record, +3 qualitative mechanism
records, Bi-2020 lane note corrected); `brands_sugar_casein_120C_pH68.yml` (+39 figure-readoff
series, header and `what_was_deliberately_NOT_extracted` amended in place);
`directional_claims_panel.yml` (+1 source, +11 claims, `hofmann1998` promoted
`access: abstract` → `full_text`); `directional_accuracy_report.md` (CURRENT STANDING
re-scored + new §W); `maillard_validation_benchmarks.md`; README, AUDIT,
VALIDATION_CONTRACT; `src/{barrier_constants,comparative_cli,curated_pathways,
directional_reliability,model_card}.py`.

**`src/model_card.collect_sulfur_anchor_status` was widened.** It inspected **one retired
file**, so what it actually measured was "the retired file still has no verified value" —
permanently true and not the question anyone is asking. It now scans the panel
(non-recursively, so the hold-out stays out) and `zero_absolute_anchors` is derived from that
scan. The guard was re-pinned to assert the **mechanism** (the flag agrees with the scan)
rather than the **answer**, which is what let it go stale in the flattering direction.

**Regenerated:** `prediction_uncertainty.{json,md}`, `matrix_binding_mode_comparison.{json,md}`,
`model_card.json` + the README card.

### (h) Gates and tests

```
holdout_guard   PASS   21 bundles · 17 maillard-path · all free_precursor · none named by any fit record
citation_gate   PASS   102 files · 925 DOI fields · 288 unique DOIs
fit_target_gate PASS   circularity + coverage accounting + constant precision
```

Targeted subsets only, per the owner's testing directive — no full-suite run.
**210 passed, 0 failed** across `test_honest_headline_guards`, `test_comparative_cli_2026_08`,
`test_wave_h_2026_08`, `test_wave_p_chemistry_2026_08`, `test_benchmarks`,
`test_free_aa_quantitative_regression`, `test_wave_s1b_ph_aw_routing_2026_08`,
`test_wave_s3_trunk_kinetics`, `test_wave_i_chemistry`, `test_wave_i_network_chemistry`,
`test_thiamine_fragmentation_benchmarks`, `test_wave_s1_additive_flux_2026_08`,
`test_wave_s4_protein_binding`.

**No PDF was committed.** Wave S5's uncommitted work was preserved (verified with
`git status` before and after; every S5 change is still staged and untouched).

### (i) [P] CARRIED FORWARD from this wave

- [P] **The step-level tables are the highest-value unexecuted evidence in the tree.** Six
  Hofmann tables constrain individual steps and none can run. Closing even one — the
  hydroxyacetaldehyde/mercapto-2-propanone → MFT row, 268.1 µg at 0.24 mol% — would give the
  repo its first single-step target. Needs species **plus** steps, not species alone.
- [P] **The binding physics is the wrong functional form for hexanal**, on the paper's own
  evidence, and hexanal is the matrix lane's most load-bearing compound. A −58.6-point
  cross-check residual is not a calibration problem.
- [P] **`K = 684.46 M⁻¹` is stranded on a missing protein molar mass**, for the second source
  running. One number — the MW the Klotz analysis assumed — would unlock both this and
  `wongprasert2024_ppi_klotz`.
- [P] **The pH sign-crossing claim is now pre-registered and unfitted.** Any future pH or
  a<sub>w</sub> term must reproduce a pentose that is monotonic and a hexose that peaks at
  pH 5, and no single family-level pH shape can do both.
- [P] **Ribose and xylose are indistinguishable to the network.** Six significant figures
  identical. The measurement separates them by 1.26–1.38×.
- [P] **The dry-heat rows need an a<sub>w</sub> basis** before they can ever be absolute rows.
  Recorded as ordinal only; do not invent one.
