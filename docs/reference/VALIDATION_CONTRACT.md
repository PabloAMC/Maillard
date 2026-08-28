# Scientific Validation Guide

If you are new to the repository, read [README.md](../../README.md) first — its architecture section and its generated model card. This document is the deeper contract and methodology reference.

## 1. Validation Contract

The project currently uses a layered benchmark contract rather than a single release number.

- The single source of truth for the benchmark contract now lives in `src/validation_contract.py`.
- **Coverage**: all measured compounds in a supported benchmark must resolve to predicted outputs.
- **Ranking**: Pearson R is only reported when at least 3 compounds match.
- **Scale**: the strict free-amino-acid gate requires full coverage plus a max measured/predicted ratio <= 1.5x.
- **Support status**: benchmarks that cannot run through the current execution path are reported as `unsupported`, not silently skipped.

The contract now separates three meanings that had started to drift together in docs and tests:

- **Directional validity**: preserves the reported ordering or trend without claiming absolute concentration fidelity.
- **Quantitative replication**: meets full-coverage plus thresholded Pearson / ratio criteria.
- **Formulation utility**: remains useful for recipe ranking or intervention comparison even when a benchmark is not strict-ready.

The benchmark-comparison policy is now explicit too:

- `free_precursor` benchmarks compare the FAST observable path, i.e. the benchmark-facing `predicted_ppb` output from the recommender.
- Cantera is currently a diagnostic reference lane for mechanism export, temporal debugging, and chemistry investigations; it is not the authoritative comparator for strict benchmark pass/fail.
- `matrix_only` benchmarks compare the dedicated intake/headspace path and therefore remain outside the strict gate and outside FAST target snapshots.

The absolute projection contract is explicit as well:

- The current projection strategy is `precursor_limited_observable_v1` in `src/recommend.py`.
- The limiting precursor pool is converted from mM to M, multiplied by a conservative severity-dependent volatile yield fraction, allocated across selected terminal outputs, and then converted to ppb on an aqueous mass-equivalent basis before matrix and headspace penalties are applied.
- Those strategy constants and assumptions are emitted in `projection_context` and in each target's `projection_metadata`.

The current benchmark summary artifact is generated with the Docker-validated helper:

```bash
./scripts/docker_maillard.sh summary
```

This writes `results/validation/benchmark_summary.md` and `results/validation/benchmark_summary.json` and is the default diagnostic surface for supported vs unsupported cases, ranking quality, and remaining scale gaps.

The benchmark index artifact is generated with:

```bash
./scripts/docker_maillard.sh index
```

This writes `results/validation/benchmark_index.md` and `results/validation/benchmark_index.json` so execution-path limits such as matrix-only benchmarks are explicit.
The index also records which engine is authoritative for each benchmark and whether Cantera is only diagnostic or not currently authoritative.

For a single graphical snapshot of reliability and current limitations, use:

```bash
./scripts/docker_maillard.sh validation-figures
```

This writes `results/validation/validation_overview.png`, `results/validation/validation_overview.md`, and `results/validation/validation_overview.json`.

For a reproducible thermodynamic-gating audit, use:

```bash
./scripts/docker_maillard.sh thermo-gating
```

This writes `results/validation/thermodynamic_gating_audit.md` and `results/validation/thermodynamic_gating_audit.json` and reports whether the gated variant materially improves benchmark error. *(2026-08-27, Wave R: both artefacts are gitignored and absent from the working tree — the command above regenerates them; this document describes what the command writes, not files you can expect to find.)* Until that audit says otherwise, thermodynamic gating remains diagnostic-only for benchmark pass/fail.
The benchmark metadata and generated summary/index artefacts now also expose the current thermodynamic-gating policy so the `auto` benchmark path resolves through an explicit contract rather than an implicit default.

For the current matrix-specific validation surfaces, use:

```bash
./scripts/docker_maillard.sh matrix-deltas
./scripts/docker_maillard.sh matrix-evidence
./scripts/docker_maillard.sh matrix-assertions
./scripts/docker_maillard.sh matrix-readiness
./scripts/docker_maillard.sh matrix-branch-deltas main
./scripts/docker_maillard.sh coverage-gaps
```

These artifacts make matrix progress and matrix limitations explicit without pretending the matrix lane is already strict-ready.

For marker-separated validation lanes, use:

```bash
./scripts/docker_maillard.sh scientific-fast
./scripts/docker_maillard.sh kinetics-validation
```

The fast lane covers FAST regression and benchmark contract checks, while the kinetics lane covers the slower Cantera reference cases, including temperature-ramp validation.
Benchmarks with fewer than 3 matched compounds now report `pass-no-ranking` when full coverage and scale pass, which keeps the summary aligned with the contract that Pearson is not meaningful for 1-2 point systems.

For reproducible target-level scientific inspection, use:

```bash
./scripts/docker_maillard.sh targets data/benchmarks/cys_ribose_140C_Hofmann1998.json
```

*(2026-08-27, Wave R: a second copy-pasteable line here ran `targets` against
`data/benchmarks/cys_glucose_150C_Farmer1999.json competing`. That file was DELETED as
fabricated — see section 3's quarantine note below — so the command could not run. It is
removed rather than repointed: there is no equivalent competing-target benchmark to
substitute, and inventing one to keep an example alive is the failure mode this document
exists to prevent.)*

This exposes the current target snapshot (`ppb`, `span`, `depth`, weighted flux) without relying on brittle inline Python or shell quoting.

For the aggregated target artifact with headspace observability metadata, use:

```bash
./scripts/docker_maillard.sh targets-report
```

This writes `results/validation/benchmark_targets.md` and `results/validation/benchmark_targets.json`, including per-target headspace class and Henry-law metadata when available.

`matrix_only` benchmarks are intentionally excluded from this target artifact. They remain executable through the summary and index artefacts, but they do not run through the free-precursor FAST target-ranking path, so exposing them as ordinary target snapshots would be scientifically misleading.

The explicit strict gate for free-amino-acid benchmarks is:

```bash
MAILLARD_STRICT_BENCHMARKS=1 ./scripts/docker_maillard.sh pytest tests/scientific/test_benchmarks.py -q
```

## 2. Current Validated Envelope

> **Provenance audit (2026-08-26).** An adversarial citation audit found that the
> sources behind three panel benchmarks are dead or resolve to unrelated papers:
> `cys_ribose_150C_Mottram1994`, `cys_glucose_150C_Farmer1999`, and
> `thiamine_cys_ribose_100C_Hofmann1996`. They are quarantined under
> `data/benchmarks/quarantined/` (excluded from the panel) pending human source
> verification — see the README there and `tasks/audit_remediation.md`. Two of
> them were previously listed below as strict-ready; that status is withdrawn.
> The same audit found the panel headline mixed literature-measured rows with
> internal synthetic comparators; coverage is now reported split by signal
> origin (see §3E).
>
> **Second provenance finding (2026-08-27, cold-start red team).** Two further panel
> benchmarks are quarantined: `spi_hvp_xylose_120C_PMC9905368` and
> `wheat_gluten_hvp_xylose_120C_PMC9905368`. **Their DOI is live and its metadata
> matches** — which is why the 225-anchor sweep passed them — but the cited paper
> (`10.1007/s10068-022-01194-w`) uses glucose/fructose at pH 7.5 for 90 min and reports
> only *relative peak areas*, and never mentions FFT or MFT. A relative-peak-area paper
> cannot be the source of an absolute ppb value. **The panel is therefore 14, not 16.**
> Read this together with §3F: the blocking citation gate checks DOI *identity* and cannot
> check DOI *content*, so this class of defect is not machine-detectable here.
>
> **Third finding, and the one that most affects how to read §3F: circularity.** A
> constant was re-derived against those two benchmarks and the same two rows were then
> scored as evidence for it — they were the only two hits in the previous "2/11 literature
> rows" headline, whose honest value was **0/11**. Coverage is now reported with fitted
> rows removed from *both* numerator and denominator, and
> `scripts/ci/fit_target_gate.py` makes undisclosed fit-then-score a build failure.

As of the current benchmark summary:

- **Strict-ready count is now 0/17 *(2026-08-28, Wave W: the panel is now 17 and strict-ready is 0/17 — three absolute Hofmann & Schieberle 1998 anchors were added and all three fail. The count of PASSES did not move.)* (2026-08-27 the panel shrank 16 -> 14 when two more benchmarks were quarantined as fabricated; 2026-08-28 it grew 14 -> 17).** The panel has no strict-ready benchmark
  at all. The sequence, all on this branch and all recorded in `tasks/audit_remediation.md`:
  6/16 → 4/16 when the projection retune cost `spi_hvp_xylose_120C_PMC9905368` and
  `wheat_gluten_hvp_xylose_120C_PMC9905368` their status, then 4/16 → 0/16 when the Wave G1
  chemistry rebuild removed the fabricated MFT shortcut and the fabricated lipid radical
  chain. `cys_ribose_140C_Hofmann1998`, the last strict-ready free-amino-acid benchmark, now
  under-predicts MFT by 5.58x and FFT by 3.26x against a 1.45x contract; the safety-marker
  benchmarks fail on the separately documented unit collisions and the SME desaturation.
  Every one of these tolerances is UNCHANGED — nothing was widened to recover the count.
  The Hofmann-only refit that followed
  (`results/validation/sulfur_barrier_refit_hofmann.md`) established that no barrier value in
  any defensible range recovers it either: the residual is a volatile-budget ALLOCATION
  deficit, not a barrier deficit.
  **2026-08-27 (Wave S2c) — the 1.45x contract named above no longer exists, and the count is
  still 0/14.** Wave S2b showed that `cys_ribose_140C_Hofmann1998`'s MFT 342 ppb and FFT
  200 ppb are a repo-internal derivation, not a measurement from `10.1021/jf9705983`: both are
  interior points of two invented, overlapping mol % bands in
  `data/benchmarks/maillard_validation_benchmarks.md` §1.3, an abstract-reconstructed table
  committed in the same commit as the benchmark JSON. The 1.45x / 0.09 dex contract was
  **~1.7x tighter than the 2.5x spread of the band its own target was interpolated from**, so
  it was **RETIRED** rather than widened; `metadata.tier` was demoted `PRIMARY → REFERENCE`,
  which removes strict-gate eligibility via `strict_gate_tiers`; and the constant that refit
  had produced (`thiol_addition_pentodiulose` 26.35) was reverted to the un-fitted class value
  28.60, its record retracted. **Strict-ready is 0/14 before and after** — this benchmark was
  *failing* its contract when it was retired (2.2086x / 0.2352 dex, and 4.3797x / 0.4041 dex
  after the revert), so retiring it removes a failure, not a pass. Retiring the contract does
  not leave the row untested: `_resolve_scale_thresholds` falls back to the global
  free-precursor defaults (1.5x / 0.10 dex), which are marginally *looser* than what was
  retired and which the row fails by more. That inheritance is stated here rather than left to
  be discovered. **The sulfur branch now has zero absolute literature anchors.**
  **CORRECTED 2026-08-28 (Wave W): THIS IS NO LONGER TRUE.** The full text of Hofmann & Schieberle 1998 (`10.1021/jf9705983`) arrived by interlibrary loan and the sulfur branch now has **three** absolute, stable-isotope-dilution literature anchors — `hofmann1998_{ribose,glucose,fructose}_cysteine_145C_20min_pH5`, the pH-5.0 aqueous rows of the paper's own Table 1 (ribose FFT 121 / MFT 198 ppb; glucose 28 / 19; fructose 32 / 25, all µg per 100 mL × 10 with the volume printed in the table footnote). The paper also confirms Wave S2b's forensic finding from the primary source: 342 and 200 ppb appear nowhere in it. **The model fails all three anchors** — 12.27×, 29.58× and 14.46× worst-ratio, mean 0.92 dex — so the branch went from *unanchored* to *anchored and measurably wrong*, which is the direction of travel this audit wanted.
- `pea_isolate_40C_PratapSingh2021` is executable through a dedicated `matrix_only` intake path, with full coverage, but it remains outside the strict release gate. Note: its marker yields are back-fitted to this same measurement, so this benchmark validates the intake plumbing, not predictive accuracy.
- `soy_isolate_40C_PratapSingh2021` is executable through the same dedicated `matrix_only` intake path, with full coverage, with the same caveat.

The headspace module now also carries a conservative matrix fallback for `pea_iso` and `soy_iso` when no explicit fat/protein fractions are supplied. That fallback is intentionally tied to the same retention estimates already used in `src/matrix_correction.py`, and those retention estimates now relax with `denaturation_state` in the same way across matrix correction, output projection, and sensory scoring.

The headspace layer now also carries the explicit observable calibration for the Pratap-Singh plant-matrix intake family: pea remains the reference yield baseline, soy-vs-pea release gaps for the tracked off-flavour markers are modeled in src/headspace.py, and the narrower pH-release correction still preserves the Pratap-Singh pH ~6 baselines while reproducing the Pouvreau acidic-vs-less-acidic pea trend.

The matrix-accessibility layer now uses explicit native/desaturated endpoints for reactive lysine and cysteine instead of implicitly relaxing all the way to free-amino-acid behavior. The code now also exposes explicit literature windows for the two isolate families it treats as quantitatively anchored: pea-isolate lysine is kept in the approximate 0.30-0.45 reactive envelope while pea cysteine remains near-zero to low-single-digit accessibility, and soy-isolate lysine/cysteine stay in a broader but still conservative 0.36-0.50 / 0.10-0.24 envelope derived from the repo's soy literature synthesis.

That calibration is now connected to a temperature/time/pH-aware denaturation heuristic in `src/matrix_correction.py`. In the formulation and optimizer paths, `denaturation_state` is no longer forced to `0.5` by default: if the user does not specify it explicitly, the engine infers an effective matrix state from the reaction conditions and exposes that `effective_denaturation_state` in the formulation result. Explicit overrides still win when the user or a benchmark needs to pin the state manually.

The user-facing explainability surface is also broader now:

- benchmark target artefacts expose proxy-vs-observable projection semantics directly (`Proxy ppb`, `Observable ppb`, `Obs/Proxy`, matrix factor, headspace factor, volatile class)
- formulation-level explainability can be generated with `./scripts/docker_maillard.sh explain-formulation <name>`
- domain-of-applicability warnings can be generated with `./scripts/docker_maillard.sh validated-envelope`

Finally, the previous bulk matrix fallback has been narrowed into a compound-class-aware overlay in the recommendation/headspace/sensory paths. Aldehydes and furans remain more strongly retained than sulfur volatiles and alcohols, which makes the default pea/soy fallback more interpretable without claiming benchmark-fitted explicit binding chemistry for every compound.

For soy, the current envelope remains conservative rather than benchmark-fitted, but it is no longer a free-floating placeholder: the repo now anchors it explicitly to the internal plant-protein literature synthesis that describes soy glycinin/β-conglycinin as compact globulins with buried reactive groups, soy as lysine-rich but sulfur-poor relative to meat-like systems, and extrusion / protein-polysaccharide conjugation as the mechanisms that partially reopen accessibility while still trapping volatiles. That is enough to justify both the relative contract `free > soy_iso > pea_iso` and the narrower calibrated window frozen in unit tests, while still keeping soy outside any strict absolute accessibility benchmark.

For `matrix_only` today, the contract is narrow and explicit:

- The benchmark must declare a non-`free` `protein_type` with a registered matrix model.
- Validation is an intake/headspace execution check, not a general precursor-resolved chemistry claim.
- It can appear as `supported` in summary/index outputs while still remaining `strict_ready = False`.
- It is deliberately omitted from `targets-report` until there is a benchmark-facing target-ranking model for matrix systems.
- The new `protein_type` headspace fallback is a conservative bridge for unspecified matrix fractions, not a replacement for explicit pea/soy composition calibration.
- Denaturation-aware retention is now internally consistent across `src/matrix_correction.py`, `src/headspace.py`, `src/recommend.py`, and `src/sensory.py`, but its numeric endpoints are still literature-estimate placeholders rather than benchmark-fitted matrix physics.
- Denaturation-aware amino-acid accessibility is now explicit too, and canonical lysine/cysteine precursor keys are recognized inside the recommendation path so matrix corrections actually apply to those species.
- That accessibility contract is now frozen end-to-end by `tests/integration/test_matrix_accessibility_recommendation.py`: sulfur-rich meaty recommendations are penalized in the order `free > soy_iso > pea_iso`, and denatured pea recovers signal without collapsing back to free-amino-acid behavior.
- The unit layer now freezes the same relative hierarchy in `tests/unit/test_matrix_correction.py`, including the stricter expectation that soy and pea concentrates stay more buried than their corresponding isolates across native, midpoint, and denatured states.
- pH-dependent plant-matrix headspace validation is now covered as a relative trend check against the Pouvreau pea-isolate family; it is intentionally narrower than the absolute Pratap-Singh intake benchmarks.

This means the framework currently has NO strict-ready benchmark (2026-08-27; see §2) and two executable matrix-headspace intake benchmarks (`pea isolate` and `soy isolate`). The matrix headspace lane should be treated as an executable intake model rather than a broadly calibrated release gate.

## 3. How We Validate

### A. Literature Benchmarks

Benchmark execution is centralized in `src/benchmark_validation.py` and reused by:

- `tests/scientific/test_benchmarks.py`
- `tests/scientific/test_free_aa_quantitative_regression.py`
- `scripts/compare_sim_to_lit.py`
- `scripts/generators/generate_benchmark_summary.py` *(2026-08-27, Wave R: path corrected; the script moved into `scripts/generators/` and this line still pointed at `scripts/`)*

### B. Sensory Fidelity

Sensory scoring remains downstream of concentration prediction. Stevens' law and masking logic are useful for ranking perceived aroma, but they do not replace benchmark validation of predicted ppb.

### C. Safety Conservatism

Safety scoring is still treated conservatively and independently from the benchmark gate so that flavor optimization does not hide formation-risk regressions.

### D. Execution Lanes

The validated execution contract now uses named Docker lanes instead of ad hoc command sequences:

- `./scripts/docker_maillard.sh core`: unit and integration correctness checks.
- `./scripts/docker_maillard.sh scientific`: benchmark summary/index generation plus scientific FAST regressions.
- `./scripts/docker_maillard.sh qm-heavy`: QM and external backend validation.
- `./scripts/docker_maillard.sh hofmann`: reproducible diagnostic trace for the calibrated Hofmann sulfur benchmark.
- `./scripts/docker_maillard.sh targets ...`: reproducible target snapshot for a specific benchmark and target family (`desirable`, `competing`, `toxic`; aliases `off_flavour` and `off-flavour`).
- `./scripts/docker_maillard.sh targets-report`: reproducible aggregate target artifact with headspace observability metadata.
- `./scripts/docker_maillard.sh validated-envelope`: reproducible validated-envelope / domain-of-applicability artifact.
- `./scripts/docker_maillard.sh explain-formulation <name>`: reproducible formulation explainability artifact for a formulation in `data/formulation_grid.yml`.

### E. External Hold-out Methodology (added 2026-08-26)

The external hold-out is the only surface that tests the frozen model on data it
has never seen. Its rules, stated explicitly:

- **Exclusion.** Hold-out bundles live in `data/benchmarks/external_validation/`
  and carry `evidence_class = external_validation_only`. The benchmark loader's
  non-recursive glob makes the directory structurally invisible to every fitting
  and panel-evaluation path, and `calibrate_from_intake` raises on the evidence
  class as a belt-and-braces guard. No hold-out value appears in any fitted
  parameter.
- **Prior tier.** Hold-out envelopes use the `uncalibrated` matrix prior tier
  (matrix_headspace ln-sigma 2.86, ~±110x at 90% CI — a figure that became true only on
  2026-08-27, when the sampler's fixed ±100x clamp was replaced by a 3-sigma-derived one;
  before that 10.7% of draws were pinned at the clamp and the realised band was ±100x)
  because their process states
  have no calibrated registry entry. This is much wider than the in-panel tier,
  so **hold-out coverage and in-panel coverage are not comparable numbers**. The
  sigma is sized from data, not judgment (2026-08-26): leave-lane-out transfer
  error over the in-panel matrix anchors — hold-out untouched — gives RMS
  ln-sigma 2.86 with a 90% interval of [1.98, 5.48] (n=6, see
  `results/validation/matrix_sigma_residual_derivation.md`). The originally
  shipped 2.0 sat at the lower edge of that interval; the tier now uses the
  residual-derived point estimate.
- **Extrapolation vs re-scoring.** Bundles whose executable conditions are
  copied from an in-panel benchmark re-score that anchor's prediction against a
  different paper's measurement — a reproducibility comparison, not an
  extrapolation test. Only bundles at genuinely new process states (HME
  extrusion, roasting) test transfer, and the honest statement is that the model
  largely fails those.
- **What the eight points ARE (added 2026-08-27).** Only **4 of the 8** are reported
  measurements. **2** are `band_geometric_midpoint` — the source reports a *range* across
  commercial lots and the scored value is `sqrt(min*max)`, a number we constructed, whose
  honest uncertainty is the 10–12x band itself. **2** are
  `derived_from_oav_and_repo_threshold` — the source's odour-activity value multiplied by
  *this repository's own* hexanal odour threshold (4.5 ppb, compilation-level, never verified
  against a primary table), so they partly encode one of our constants and **move if we
  correct it**. Every row carries its `value_provenance` and the report renders the split.
  Until 2026-08-27 these were written at full float precision beside a fabricated
  `measurement_date` of "<publication year>-01-01", which is now `not_applicable`.
- **Current numbers (2026-08-27 Wave R regeneration). The headline is 1/5.**
  Genuine-extrapolation coverage at the **pre-widening** prior (ln-sigma 2.0) is
  **1 of 5**, and since Wave O it is 1/5 under the shipped ln-sigma 2.86 as well (it was
  2/5 between Wave M and Wave O; that difference was the width of the interval, not the
  accuracy of the model). Over all eight rows: **3/8 at ln-sigma 2.0 and 3/8 at 2.86** — the
  two priors now agree, because Wave R's Liu correction pushed the last point that separated
  them out of both intervals. Of the covered rows **2 of 3 are re-scoring bundles**
  (`bi_2020_raw_pea`, `liu_2023_ppi_offnote_baseline`) that test nothing. Median fold error
  **93.68x** (median |log10| 1.972 dex); worst miss **2474x** on roasted pea. **Read the fold
  error, not the coverage** — it is the only figure here no choice of prior can move.
  **Wave R (2026-08-27) made this number WORSE by replacing two reference values that were in
  no source.** The `liu_2023_ppi_offnote_baseline` bundle carried hexanal 15–180 ppb and
  nonanal 5–50 ppb attributed to "Liu, Y. (2023 thesis)". The thesis (Yaozheng Liu, *Flavor
  Chemistry of Pea Proteins*, NC State 2021; published as Liu, Cadwallader & Drake 2023,
  Food Chem. 406:134998, `10.1016/j.foodchem.2022.134998`) was retrieved and read in full:
  Table 2.7 reports hexanal **2445–52454 µg/L** and nonanal **0.188–3.42 µg/L** across nine
  quantified commercial pea proteins, and no table in the document contains the repo's bands
  or the OAV pairs attached to them. Correcting the measured values — **no prediction moved,
  nothing fitted** — gives hexanal 19.50x → **11.17x** and nonanal 4.78x → **94.22x**; median
  42.62x → 93.68x, coverage 4/8 → **3/8**, `max_fold_error` and the pre-widening 1/5
  unchanged. The nonanal row is now the sharpest lipid-lane over-prediction the repository
  has against a directly-quantified reference (75.5 ppb predicted, band top 3.42 ppb); Wave
  P's oleate-substrate fix is the partial mitigation already landed and took the same point
  from 214x to 94x. Liu's standard curves were built in **deionized water, not in the protein
  matrix**, so protein binding is uncorrected and her values are lower bounds — the gap is if
  anything understated. Two further anchors on the same source were **retired, not repaired**:
  (E,E)-2,4-heptadienal is absent from the thesis entirely, and the methoxypyrazine anchor
  named IBMP, which the thesis neither identifies nor quantifies (its methoxypyrazine is IPMP
  at 6.126–57.0 µg/L, 713x the retired 0.08 ppb ceiling). Neither was ever in the scored
  bundle. The hold-out stayed out of every fit: this was a reference correction, not a
  calibration.
  **Wave O (2026-08-27) made this number WORSE by making a constant more correct.** The
  ambient hexanal observability factors were refitted from the 260 / 380 ppb transcription
  errors onto the paper's verified 1138.00 / 1621.71 ppb (one shared scale, 4.317249x;
  record `results/validation/matrix_observability_refit_pratap_singh.json`). Median 15.31x
  → 42.62x, coverage 5/8 → 4/8. Bi 2020 raw pea improved 5.37x → 1.24x; Liu 2023 went
  4.52x → 19.50x and Li 2026 hexanal 21.58x → 93.15x. The reason is a contradiction inside
  the hold-out itself: Bi reports 1260 ppb and Liu's band midpoint is 51.96 ppb for
  nominally the same system, a 24x spread, and the erroneous 260 ppb sat almost exactly at
  their geometric mean (255.9). No observability factor satisfies both. `max_fold_error`
  and the pre-widening 1/5 did not move.
  *(Superseded 2026-08-27 by Wave R: the 24x spread was not a literature contradiction. The
  51.96 ppb midpoint came from a band that appears in no source; Liu's real Table 2.7 range is
  2445–52454 µg/L, which puts the verified 1138 ppb anchor just under her lowest lot. The Wave
  O refit's direction is vindicated by the corrected target — this same point reads 11.17x
  against the real number.)*
  **Wave P (2026-08-27) improved two of the eight points with nothing fitted, and moved the
  headline by nothing.** Nonanal is the C9 fragment of the OLEATE double bond, not a
  linoleate product (Miyazaki 2023, 10.1093/bbb/zbac189; Hung, Katrib & Martin 2005,
  10.1021/jp0500900), and `LipidProfile.oleic_acid_pct` had been dead code. Correcting the
  substrate moved `li_2026_spi_wg_hme_control` nonanal 272.63x → **118.31x** and
  `liu_2023_ppi_offnote_baseline` nonanal 10.86x → **4.78x**, each by exactly its matrix's
  oleic/linoleic ratio; the other six points are byte-identical. (Both Liu folds here are
  against the reference value Wave R later retired; against Table 2.7 the same fix reads
  214x → 94.22x.) Median fold error, coverage
  hits (4/8, now 3/8 after Wave R), `max_fold_error` (2474x) and the pre-widening 1/5 were
  **all unchanged by this fix** — the
  median sits between two points that did not move, and the two that did were outside the
  interval before and after. This is the clearest illustration in the repository of why the
  per-point table has to be read alongside the headline. Both nonanal points remain
  over-predicted, and the model still treats oleate as being as oxidisable as linoleate.
  **Wave S1 (2026-08-27) left all eight points BIT-IDENTICAL, and the reason is a finding in
  its own right.** Two structural fixes landed in `src/recommend.py` -- the flux propagator
  became additive over parallel channels, and the compound-specific matrix calibration
  registry became reachable on the `matrix_precursor_augmented` lane -- and together they
  moved 26 of the 42 scored in-panel rows. They moved ZERO of the eight hold-out points.
  Every hold-out bundle executes the `matrix_only` path, which passes compound NAMES to the
  registry (so the lookup repair never applied to it) and bypasses `predict_from_steps`
  entirely (so the propagator change never applied to it). **The external hold-out therefore
  exercises the lipid-oxidation and observability lane and says nothing whatsoever about the
  Maillard network propagator.** Median 93.68x, coverage 3/8, `max_fold_error` 2474x and the
  pre-widening 1/5 are all unchanged, and that invariance is evidence about the hold-out's
  coverage rather than evidence about the model.
  **Wave S1b (2026-08-27) left them BIT-IDENTICAL AGAIN, for the same structural reason.**
  Three pH / water-activity ROUTING defects were repaired in `src/conditions.py` -- the
  enolisation route-selection term had never been called on the prediction path, the pyrazine
  ionisation branch keyed on a substring matching none of the 29 emitted families, and the
  water-activity correction reached 3 of those 29 and missed the furan track. Those changes
  moved 4 of the 14 benchmark rows (all four in the wrong direction), moved the directional
  panel's pH bucket from 2/7 to 4/7, and moved the MC panel's honest literature coverage from
  1/3 to **0/3**. They moved ZERO of the eight hold-out points, because `get_rate_constant`
  is only reached through `predict_from_steps`, which `matrix_only` bypasses. **Two
  consecutive waves of prediction-path changes have now been invisible to the external
  hold-out.** Until a hold-out bundle exercises the Maillard network, no external evidence
  bears on the propagator, the barriers, the network topology, or the pH/aw physics.
  **Wave S2c (2026-08-27) left them BIT-IDENTICAL A THIRD TIME, and this one is the sharpest
  illustration of the gap.** Wave S2c reverted a shipped BARRIER —
  `thiol_addition_pentodiulose` 26.35 → 28.60 kcal/mol — after Wave S2b established that the
  benchmark it had been fitted against, `cys_ribose_140C_Hofmann1998`, carries values this
  repository derived rather than measured (interior points of two invented mol % bands in
  `data/benchmarks/maillard_validation_benchmarks.md` §1.3, committed in the same commit as
  the benchmark JSON). That revert moved the flagship MFT row by ~2.0x, moved the
  pentose ≫ hexose ordering headline from 18.27x to 8.26x, and re-scored two panel
  benchmarks. It moved **zero** of the eight hold-out points, for the same structural reason:
  `matrix_only` never reaches `predict_from_steps`, so no `FAST_BARRIERS` entry is on its
  path. **Three consecutive waves — a propagator change, a routing repair, and a barrier
  revert — have each been completely invisible to the external hold-out.** That matters more
  after Wave S2c than before it: the hold-out cannot corroborate or refute a single
  sulfur-branch constant, and the sulfur branch now has **zero absolute literature anchors**
  in-panel either. There is currently no evidence of any kind bearing on those constants.
  **CORRECTED 2026-08-28 (Wave W): the second half of that is no longer true.** The panel
  now carries three absolute stable-isotope-dilution sulfur anchors from Hofmann &
  Schieberle 1998, and five more rows from the same paper sit in the Maillard-path
  hold-out. There is now evidence bearing on the sulfur constants, and it is unflattering:
  the shipped barriers miss the three in-panel anchors by 12.27x, 29.58x and 14.46x. What
  is STILL true is the first half — the *matrix-only* hold-out remains blind to every
  sulfur constant, because it never reaches `predict_from_steps`.
  **The 0/5 → 1/5 and the 32.79x → 15.31x are a REFERENCE correction, not a model change
  (Wave K/M, 2026-08-27).** Two of the four `li_2026_spi_wg_hme_control` points had been
  transcribed from adjacent table rows: 2-pentylfuran 221.5 was the paper's *Maltol* row
  (true 5625.80 ppb) and nonanal 29.42 was its *Decanal* row (true 72.66 ppb), verified
  against Europe PMC `fullTextXML` (PMC12984281). The 2-pentylfuran point moves 49.8x over
  → 1.96x over and nonanal 673x → 273x with byte-identical predictions. The extreme
  process-state misses (roasted pea 2474x, HME 1-hexanol 1117x) are untouched. Historical:
  33.84x before the sigma was raised, 32.79x after the Wave I regeneration.

### E-bis. Curating a NEW lipid-oxidation hold-out anchor

*Folded in 2026-08-28 (Wave S5) from `docs/guides/CURATING_LIPID_OXIDATION_ANCHOR.md`, which was deleted. It was current and correct but orphaned — nothing linked it except the remediation ledger — and it is methodology, which is what this document is for. Its `## Related` tail is dropped: one of its three links pointed at `docs/architecture.md`, which this wave also deleted. The live diagnostic is still `scripts/diagnose_lipid_bias.py`.*

> **Gate (hard rule): do not invent numbers.** This checklist describes *how* to add a
> real, literature-sourced lipid-oxidation anchor to the frozen external hold-out. Do not
> add an entry until a qualifying paper is in hand. The 4 existing bundles in
> `data/benchmarks/external_validation/` remain the frozen test set.

#### Why this matters

S27 established that the matrix lipid-aldehyde over-prediction on the external hold-out is a
**per-(matrix, process_state) calibration gap**, not a kinetic-shape problem
(see [src/lipid_oxidation.py](../../src/lipid_oxidation.py) and
[data/lit/lipid_oxidation_calibration.json](../../data/lit/lipid_oxidation_calibration.json)).
Workstream B made the model *honestly uncertain* on uncalibrated process-states (external
hold-out 0/8 → 3/8 inside 90% CI; 5/8 as of 2026-08-27, after the uncalibrated matrix
ln-sigma was raised 2.0 → 2.86 on residual evidence — a wider interval, not a better
prediction: median hold-out fold error is ~33× either way). The extreme-processing cases
(HME extrusion, roasting) remain genuine misses. **Closing them requires real measured anchors at those process-states**,
not a wider prior or a tuned cap. That is the curation task below.

> **Update 2026-08-27 (Wave G1/H chemistry rebuild).** The network's lipid radical chain was
> rebuilt: `radical_propagation_o2` had been matching *any* sp² carbon, so 61 of 103
> peroxidation steps were fabricated, radical flags were being lost in the SMIRKS ordering,
> and the β-scission pattern required an sp³ β-carbon — which made **hexanal structurally
> unreachable in the network**, masked in production only by the fixed branching ratio. All
> of that is fixed and hexanal is now reachable. Two things follow for this checklist.
> (1) The chain initiates only from a **hydroperoxide**: an unoxidised fatty acid plus O₂
> enumerates to zero steps, so the anchor still seeds the lane and the calibration registry
> is still where the quantitative content lives. (2) The hold-out numbers above are
> unchanged by the rebuild (median fold error 32.79×, worst 2474× on roasted pea), because
> the hold-out lanes run through the lipid-oxidation anchor rather than the network. The
> curation task is unchanged.

#### What qualifies as an anchor (intake criteria)

A candidate paper must report **all** of the following for a well-defined plant-protein matrix:

1. **Quantitative** headspace concentrations (ppb or convertible) for at least one lipid-oxidation
   marker: hexanal, nonanal, 2-pentylfuran, 1-hexanol, or (E,E)-2,4-decadienal.
2. **SIDA-grade or internal-standard-calibrated** quantification (HS-SPME-GC-MS/SIM acceptable;
   semi-quant peak areas are **not** acceptable for a hold-out anchor).
3. A **defined matrix** (protein type + isolate/concentrate) and **process-state** parameters:
   temperature, time, and (ideally) water activity / moisture regime / SME for extrusion.
4. A process-state that **fills a calibration gap** — i.e. roasting, HME extrusion, or another
   high-severity regime not already pinned by the calibration registry. (Ambient slurry is
   already represented.)
5. Enough method detail to set `analytical_context` (method, quantification_mode, replicates).

If any of 1–3 is missing, the paper is **not** hold-out-grade — log it in the backlog instead.

#### Steps to land an anchor (once a qualifying paper is in hand)

1. **Extract** the measured marker concentrations and the full process-state, in the paper's own
   units. Record the DOI. Convert to ppb explicitly; show the conversion.
2. **Add a bundle spec** to `_HOLDOUT_BUNDLE_SPECS` in
   [src/external_validation.py](../../src/external_validation.py), mirroring the existing entries
   (`bundle_id`, `anchor_ids`, `matrix_context`, `protein_type`, `process_state`, `conditions`,
   `precursors`, `benchmark_alignment`, `analytical_context`). Keep
   `evidence_class = external_validation_only` so it stays out of calibration.
3. **Generate** the benchmark JSON under `data/benchmarks/external_validation/` via the existing
   payload generator (`scripts/generators/generate_external_validation_payloads.py`).
4. **Regenerate** the report: `./scripts/docker_maillard.sh run "python scripts/generators/generate_external_validation_report.py"`.
   The new anchor will appear with its honest fold-error and CI status — **report it as-is**.
5. **Do not** retune `max_conversion_fraction`, the observable sigmas, or any prior to make the new
   anchor pass. Anchors are a test set, not a fit target. If the new point is a large miss, that is
   a true signal that the process-state needs *dedicated calibration data*, which is a separate,
   in-panel ingestion task (`./scripts/docker_maillard.sh ingest`).

#### If the paper is good but not hold-out-grade

Log it in the literature backlog (`results/validation/literature_backlog.*`) with the reason it
was rejected (e.g. semi-quant only, undefined matrix). It may still inform the in-panel
calibration once a quantitative companion is found.

### F. Reading The Coverage Split

`results/validation/prediction_uncertainty.md` now reports coverage split by
signal origin:

- **External literature** rows are validation evidence.
- **Fitted rows** (added 2026-08-27) are rows whose *own constants were solved from them*,
  with enough freedom to reproduce them row by row. Agreement there is algebraic recovery,
  not prediction, so they are removed from the literature numerator **and** denominator and
  reported separately. **Read their outcomes rather than skipping them**: a row the model
  still fails *after being fitted to it* is a stronger negative result than any coverage
  number. `src/fit_target_index.py` classifies each fit record by leverage — a single global
  constant fitted across two dozen rows does **not** trigger exclusion, because excluding
  those rows would delete genuine failures instead of exposing them.
- **Internal synthetic** rows compare the model against its own frozen output
  (drift detection); they carry zero validation weight and exist so refactors
  cannot silently change predictions.
- **Not evaluable** rows have degenerate (near-zero-width) envelopes: the Monte
  Carlo perturbs nothing on their path, so pass/fail is meaningless for them.
- Coverage is only interpretable next to median CI width: a 90% interval
  spanning several orders of magnitude makes coverage cheap. Report both,
  always.

**Current in-panel numbers (2026-08-27 Wave I regeneration, n=200, seed 0).** **14**
benchmarks (two quarantined as fabricated since the previous revision); 11 in the Monte
Carlo panel, 35 matched rows; aggregate 90% CI coverage 29/35 (82.9%) — **do not quote that
number**, it pools three populations that support different claims:

| Population | Inside 90% CI | Not evaluable | Median CI width | Evidence? |
| --- | ---: | ---: | ---: | --- |
| External literature | **1/3 (33.3%)** | 4 | 0.86 dex | **yes — only this row** |
| Fitted rows | 2/2 (100%) | 0 | 2.28 dex | no — algebraic recovery |
| Internal synthetic | 18/18 (100%) | 8 | 3.65 dex | no — reproducibility harness |

The literature denominator fell from 11 to 3 because the two rows that used to be its *only*
hits were fitted rows, and two more benchmarks were quarantined. Benchmark-level: **0/6
predictive** benchmarks are without blocking gaps (3/4 fit-recovery and 4/4 synthetic are —
i.e. every "pass" in the panel sits in a non-evidence bucket). Strict-ready: **0/17** *(2026-08-28, Wave W: the panel is now 17 and strict-ready is 0/17 — three absolute Hofmann & Schieberle 1998 anchors were added and all three fail. The count of PASSES did not move.)* (see §2).

Read every coverage figure with its width column. Note also that the Monte-Carlo intervals
here are **wider than any previously published**, for a reason that is a defect and not an
improvement: until 2026-08-27, 10 of the 14 barrier-family priors resolved to keys the engine
never emits, so ~70% of the barrier channel was inert and every interval was narrower than the
priors claimed. Coverage numbers from before that date were computed against intervals the
model was not actually entitled to.

## 4. How To Read The Validation Figures

The repository now relies on a small set of figures that each have a distinct purpose.

- `validation_overview.png` is the first-pass trust surface. It is now a single-panel parity figure for the authoritative benchmark points, with formatted study references instead of raw benchmark ids, so it can be embedded directly in the README without multi-panel cropping.
- `validated_envelope.png` is the boundary figure. It shows how much of the benchmark set is strict-ready, how much is only directional, and which caveats still define the current edge of reliable use.
- `*_comparison.png` files are benchmark cards. They show parity, absolute yields, and the benchmark summary for one literature system at a time.

Read them in that order. The overview tells you whether the repo has a real quantitative proof surface. The envelope tells you where that proof stops. The benchmark card tells you why a single literature case passed, failed, or remained directional.

Correlation is informative but not sufficient on its own. The contract should always be interpreted through four signals together:

- coverage
- Pearson $R$ only when at least three compounds match
- MAE in ppb
- max ratio as the release-facing scale tolerance

This is why a benchmark may be scientifically acceptable as `pass-no-ranking` while still lacking a meaningful correlation coefficient.

## 5. Blind Spots That Still Matter

- **Citation provenance**: literature anchors were ingested by LLM-assisted deep-research campaigns and are not yet human-verified. An automated sweep (2026-08-26) found ~20% of registry DOIs unresolvable, plus a harder-to-detect class of live DOIs pointing at the wrong paper. Kinetic parameters and benchmark values should be treated as provisional pending manual source verification; `data/lit/computational_priors.json` (`dft_kinetic_priors`) is the provenance template the rest of the registries should be held to.
- **Matrix-only systems**: the first plant-isolate benchmark path is now executable and the default matrix-state handling is less manual, but broader matrix calibration beyond the current pea/soy intake family is still missing.
- **Headspace translation**: the remaining free-amino-acid scale gaps are now dominated by how FAST activity is translated into observed concentration/headspace, not by a single missing sulfur barrier.
- **Peptide accessibility**: intact protein systems are still outside the validated free-precursor envelope.
- **User-facing explainability**: the engine now computes a more physical matrix state internally, but that matrix metadata is not yet surfaced broadly in benchmark summary/report artifacts.
- **Default CLI integration**: explainability and validated-envelope artifacts exist, but the default user commands do not yet surface those warnings inline during every formulation run.

## 6. Recommended Verification Workflow

1. For macOS/OrbStack or Docker Desktop, bring up the validated Linux environment with `./scripts/docker_maillard.sh up` and create or refresh the env with `./scripts/docker_maillard.sh bootstrap`.
2. Run `./scripts/docker_maillard.sh summary` to inspect the current validated envelope.
3. Run `./scripts/docker_maillard.sh index` to inspect benchmark tiering and execution-path scope, especially for matrix-only benchmarks.
4. Run `./scripts/docker_maillard.sh pytest tests/scientific/test_free_aa_quantitative_regression.py -q` to guard the currently calibrated free-amino-acid ratios.
5. Run `./scripts/docker_maillard.sh pytest tests/` for the full Docker-validated suite.
6. Run `MAILLARD_STRICT_BENCHMARKS=1 ./scripts/docker_maillard.sh pytest tests/scientific/test_benchmarks.py -q` when you want an honest go/no-go signal for the strict free-AA gate.
7. Run `./scripts/docker_maillard.sh targets data/benchmarks/<benchmark>.json` when you need the current target snapshot for branch-level analysis.
8. Run `./scripts/docker_maillard.sh targets-report` when you need the aggregate target artifact regenerated before review or comparison.
9. Treat `matrix_only` benchmarks as executable intake checks unless they are explicitly promoted into the strict gate.

## 7. Expected Skips In The Docker Lane

- `tests/benchmarks/` **no longer exists.** *(Corrected 2026-08-27, Wave R.)* This line used to describe it in the present tense as "intentionally skip-heavy today … Phase 3 placeholders and HPC-oriented literature checks". The directory and its `_lane_policy.py` gate were deleted on 2026-08-27 (Wave J2); nothing skips there because nothing runs there. The loader `src/authority_benchmark_data.py` and its tracked fixtures (`data/qm/phase33_barrier_benchmarks.json`, `data/qm/phase35_double_hybrid_benchmarks.json`) survive, but the loader's only remaining consumer is its own parse test, so no value from that lane reaches the model, a calibration or a headline. Reviving or retiring the Phase 3 authority lane is an open decision, not a skip.
- Capability-gated QM tests should skip only when the backend is genuinely unavailable or unusable in the active Docker environment.
- A path-based skip for a binary that is actually present in `PATH` is a test bug, not a valid environment gate.

## 8. Named Docker Lanes

- `./scripts/docker_maillard.sh core`: unit and integration correctness gate.
- `./scripts/docker_maillard.sh scientific`: benchmark summary/index plus scientific regression lane.
- `./scripts/docker_maillard.sh qm-heavy`: QM and external backend lane.
- `./scripts/docker_maillard.sh hofmann`: diagnostic trace for the calibrated Hofmann sulfur benchmark.
- `./scripts/docker_maillard.sh targets ...`: benchmark target snapshot lane for scientific inspection.
- `./scripts/docker_maillard.sh targets-report`: aggregate target artifact lane for scientific inspection and review.

These lanes keep infrastructure validation separate from unresolved scientific scope gaps such as matrix-only precursor systems.
