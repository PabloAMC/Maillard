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
