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
