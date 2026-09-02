# The Audit

*2026-08-26/27. This document is the narrative record of a two-day adversarial audit and
remediation of this repository. The machine-readable ledger is
[tasks/audit_remediation.md](tasks/audit_remediation.md); the evidence artifacts it cites
are tracked under [results/validation/](results/validation/).*

> **Status, 2026-09-03 (retirement step B5).** The engine this audit examined — the SMIRKS
> rule-enumeration screening lane with its fitted volatile budget, its validation harness
> (`src/benchmark_validation.py`), its Monte-Carlo sampler and its matrix/headspace layer — has
> been **deleted**. The kinetic core (`src/kinetic_core/`) is the only engine, scored on its own
> (`results/validation/core_panel_scores.md`, `core_prediction_uncertainty.md`; pinned by
> `tests/scientific/test_core_headline_guards.py`). Every number quoted below was produced by the
> retired lane and is kept as the record of what that lane claimed; its artifacts are archived
> under `results/legacy_lane/` and its README under `docs/history/`. The findings about
> *process* — circular validation, fabricated citations, fit-then-score — are the reason the
> core carries a pre-registration per fit wave, a fit-target ledger and five CI gates.

## Why this document exists

This codebase was built with heavy LLM assistance. In August 2026 the owner commissioned
an adversarial audit before showing the project to external scientists. The audit found
serious problems. Rather than quietly fixing them, this repository publishes the findings,
the fixes, and the numbers that got *worse* when fabricated support was removed — because
the project's only defensible claim is calibrated honesty, and that claim is only credible
if the audit trail is public.

## What the audit found

Four independent adversarial reviews (calibration-leakage forensics, citation
verification, dimensional analysis, an end-to-end prediction trace), followed by a
mechanism-level chemistry review and a module-by-module bug hunt. The major findings, all
since remediated:

1. **The validation headline was substantially circular.** 26 of 48 "matched literature
   data points" were the model's own frozen output re-ingested as measurements — including
   files labeled as internal experiments, with measurement dates and instrument
   uncertainties, whose values were bit-identical to model output at 17 significant
   figures. The reference calibration lane reproduced its own anchor by construction.
2. **The literature anchor layer was heavily contaminated.** A full sweep of 225 anchor
   DOIs found 61% defective: 20% dead, 20% resolving to entirely unrelated papers (froth
   flotation, strawberry harvests), 21% metadata-mismatched. Four panel benchmarks cited
   sources that do not exist; their values were invented (one reported MFT = 3.14 ppb —
   π). The benchmarks with the most suspect provenance carried the *tightest* tolerance
   contracts — the signature of values fitted to model output.
3. **The flagship chemistry was fabricated in places.** The reaction network produced its
   headline compound (2-methyl-3-furanthiol) through a one-step shortcut that skips the
   sugar-fragmentation mechanism entirely; a fictitious hydrogen economy manufactured H₂ from
   water; the lipid radical chain was non-functional (93% of its steps were artifacts of
   one over-permissive SMARTS pattern, while hexanal — the off-note the model exists to
   predict — was structurally unreachable); eight reaction families, including the
   acrylamide and CML/CEL safety lanes, were silently disabled by a default-barrier
   fallthrough.
4. **The safety panel's "validated" status was a unit-collision artifact** (mg/kg values
   stored in ppb fields — a 1000× error yielding near-perfect apparent agreement), the
   sulfur barrier constants had been jointly fitted to the fabricated benchmarks (commit
   forensics), a shipped results database silently overrode curated constants depending on
   the working directory, and the projection layer applied Boltzmann selectivity twice at
   an effective temperature of T/2.19.

## What was done

Roughly 30 agent-executed workstreams over two days, all recorded in the ledger:

- **Provenance**: every defective anchor was repaired with a CrossRef-verified source or
  explicitly labeled `no_verifiable_source`; fabricated benchmarks were quarantined and
  then deleted after source-recovery research confirmed no source exists; one benchmark
  was rebuilt from a verified 1994 source — which the model honestly fails at ~773×. A CI
  citation gate ([scripts/ci/citation_gate.py](scripts/ci/citation_gate.py)) now blocks
  contamination regrowth, at zero waivers.
- **Chemistry**: a mechanistic MFT route replaced the one-step shortcut (via norfuraneol in
  Round 1 — see Round 3 for why *that* route was retired too, and note that the parenthetical
  this bullet originally carried, "the pentose≫hexose selectivity is now structural, not
  fitted", was **wrong when written**: only ~1.1× of it is structural); the fabricated radical chain was removed
  (5500 → 369 steps in the full network) and hexanal made genuinely reachable; wrong
  structures (the flagship disulfide, both furanones) were corrected against InChIKeys;
  every emitted family now has an explicit, honestly-tiered barrier.
- **Physics**: the double Boltzmann was removed and the volatile budget re-derived as an
  Arrhenius conversion extent (eliminating a spurious furfural optimum at ~145 °C that
  was an artifact of a saturating sigmoid); the Henry's-law table was rebuilt against the
  Sander 2023 compilation; the QRRHO entropy correction was fixed (it had the sign of its
  own purpose inverted) with version-guarding against stale database rows.
- **Honest reporting**: coverage is split by signal origin (literature vs synthetic
  drift-detection rows); zero-width envelopes are excluded as not-evaluable; the hold-out
  prior is sized from a leave-lane-out residual derivation
  ([results/validation/matrix_sigma_residual_derivation.md](results/validation/matrix_sigma_residual_derivation.md))
  instead of judgment; the hold-out methodology is documented in
  [VALIDATION_CONTRACT.md](docs/reference/VALIDATION_CONTRACT.md) §3E.
- **Everything else**: the safety layer got a units contract and a real [0,1] risk score;
  the optimizer now reports the formulation it actually optimized; the front-door
  commands do what the README says; runs are bit-deterministic; the dead code island was
  deleted; CI runs the suite and the integrity gates.

## Where the numbers stand (2026-08-27)

These are worse than the numbers this repository advertised a week ago. They are real.

| Surface | Value |
| --- | --- |
| Panel | **23 benchmarks** *(2026-08-28, Wave X: 17 -> 23; six STEP-LEVEL isotope-dilution rows added from Tables 3, 4, 8 and 10 of the same paper. 2026-08-28, Wave W: 14 -> 17. 2 more remain quarantined from Round 2, and **nine further verified rows are deliberately OFF the panel** in `data/benchmarks/step_level_unreachable/` because the network cannot execute them — scoring a structural zero would cost nothing, since the log-error aggregate skips non-positive predictions)*; the MC uncertainty panel is regenerated with them |
| Literature rows inside 90% CI | **4/13 evaluable** on a median CI width of **1.44 dex** — a factor of ~28 end to end. *(2026-08-28, Wave X: 4/9 at 2.63 dex -> 4/13 at 1.44 dex. Read this the opposite way round from Wave W's move: the interval got NARROWER by a factor of ~15 end to end AND the coverage rate FELL, 44% -> 31%. The four new step-level literature rows carry envelopes of 0.15–0.4 dex and every one of them MISSES. The model is not merely wrong on a single step; it is CONFIDENTLY wrong.)* *(2026-08-28, Wave W: 0/3 at 0.75 dex -> 4/9 at 2.63 dex. The rate rose because six new external-literature rows arrived, four of which sit inside intervals 2.9-3.7 dex wide; an interval spanning three orders of magnitude covers almost anything. NO prior, threshold or interval was widened in Wave W — the population changed, not the model. The two rows that MISS, fructose FFT and MFT, are the informative half.)* Fitted rows are removed from numerator and denominator alike |
| Fitted rows (constants back-solved from the benchmark) | 2/2 — **not evidence**; both would previously have been counted as literature hits |
| Internal-synthetic rows (reproducibility only) | 18/18 — carries zero validation weight |
| Benchmarks without blocking gaps | **0/14 predictive** (+ **0/5** fit-recovery, 4/4 synthetic; the **4/23** aggregate is all non-evidence). *(2026-08-28, Wave X: predictive 0/9 -> 0/14 and fit-recovery 0/4 -> 0/5. Five of the six new step-level rows are predictive and all five FAIL; the sixth is a declared fit target and moves to fit-recovery automatically. Moving a row into fit-recovery normally flatters the headline by shrinking the predictive denominator — it does not here: the denominator grew by five, the numerator stayed at zero, and the excluded row is itself a 2.3x miss.)* *(2026-08-28, Wave W: predictive 0/6 -> 0/9. All three new Hofmann anchors are predictive and all three FAIL — 12.27x, 29.58x, 14.46x worst-ratio. The headline denominator grew by 50% and the numerator stayed at zero.)* Fit-recovery fell 3/4 -> 1/4 in Round 3 when the two Pratap-Singh benchmarks were corrected against the paper and stopped recovering. Wave O then refitted their constants onto the verified values and they recover again — at `pass-no-ranking`, which is not `pass`, so none of those counts moved. A refit changes the size of a recovery, never its evidential status. |
| External hold-out | **1/5 genuine extrapolations** at the pre-widening prior — and 1/5 under the wider one too since Wave O; median fold error **93.68×**, worst **2474×**, coverage **3/8** at the shipped sigma. Only 4 of the 8 points are measurements at all. The median was 32.79× until Round 3 corrected two wrong table rows in the Li 2026 bundle (→ 15.31×, predictions unmoved), rose to 42.62× when Wave O refitted the ambient hexanal observability onto the paper's verified anchor, and rose again to 93.68× when Wave R found the Liu hold-out's two reference values matched nothing in their source and replaced them with the thesis's own Table 2.7. **Every rise came from correcting data toward the literature, and predictions moved in none of them**; see Round 3, Wave O and Wave R. |
| Strict-ready benchmarks | **0/23** *(2026-08-28, Wave X: 0/17 -> 0/23. The six step-level rows are not strict-ready either. They are predicted BETTER than the end-to-end rows — 0.518 dex against 0.952 — and still nowhere near the +/-10% replicate contract their own source's footnote sets)* *(0/14 before Wave W; the three Wave W literature anchors are not strict-ready either — they fail their contracts by 12-30x)* |
| Kinetic core, scored on its own *(2026-09-03, retirement B3; the lane that replaces the one above)* | 40 bundles / 49 rows: **8/49 within 3x**, out-of-sample **4/40**; envelope **11/44** evaluable literature rows inside the 90% CI (B8 Laplace covariance samples the sulfur lane; 5 rows not evaluable); **strict-ready 1/40** (`thiamine_cys_glucose_120C_Bolton1994`, under its own 3x contract). The sulfur fit read 8 panel rows including the xylose pH-5 hold-out. `results/validation/core_panel_scores.md`, `core_prediction_uncertainty.md` |
| Reaction-chemistry lanes with generative templates | 5/16 (derived by enumeration, test-pinned) |
| Test suite | **1274 passed, 1 skipped, 2 xfailed, 0 failed** (2026-08-27, after Wave S1; 1265 after Wave P, 1242 after Round 3). The dvipng failure carried as "environmental" for the whole audit is **gone**: `--report` now degrades to mathtext instead of raising |
| Citation gate | **0 dead DOIs as of the 2026-08-26 and 2026-08-27 sweeps; 0 waivers.** The blocking gate is **structural and offline** — DOI grammar, confabulation signatures, status coherence, repair-record completeness. It cannot detect a live DOI that resolves to the wrong paper, which is exactly how the two PMC9905368 benchmarks survived every previous check. Liveness is a **scheduled weekly sample**, advisory, not part of the required set. |

The coverage numbers *fell* during remediation — twice — and each fall is documented at
the point it happened (see the README's calibration section and the ledger). That is the
audit working as intended: each drop is the measured size of a fabricated support that was
removed. What survives at quantitative fidelity is narrow; what survives robustly is
**ordering** — pentose ≫ hexose sulfur chemistry (**8.26×** against a 3× floor; published as
15.8×, then corrected to **8.98**× on 2026-08-27, then to 3.39× the same day when the MFT
route was corrected on isotope evidence, then to 6.15× when Wave P refitted the sulfur-addition
barrier, then to 7.78× when Wave S1 made the flux propagator additive, then to 18.27× when
Wave S1b reconnected the water-activity physics, then **down to 8.26×** when Wave S2c reverted
that Wave P refit as circular; and of that, ~4.27× is
structural — zeroing the now-1.05 kcal/mol barrier gap between the hexose and pentose
sulfur-addition steps collapses it, so ~1.9× is still carried by a barrier difference rather
than by mechanism. The structural share in absolute terms has gone 1.13× → 2.31× → 3.14× →
4.27× → 4.27×, but neither the S1b rise nor the S2c fall is a sugar-chemistry finding: under
Wave S1b, at aw 0.98 the hexose limb's three-water shortcut is penalised harder than the
pentose limb's two steps, so **both** absolute numbers fell (ribose 824.7 → 374.0 ppb, glucose
106.0 → 20.5 ppb) and the ratio rose because the denominator fell faster; under Wave S2c only
the pentose limb moved (ribose 374.0 → 169.1 ppb, **glucose unchanged at 20.5**) because only
it runs the reverted barrier. **The claim shrank and its evidence improved at the same time:**
the share carried by mechanism rather than by a barrier gap went from 23% to 52%, and no part
of it now traces to a value this repository invented. **Half the claim is still a
gap between two barriers**), matrix directionality — plus an experiment-prioritization machinery whose gap map is now honest.

## Round 2: the cold-start red team (2026-08-27)

The remediation above ended by commissioning two **zero-context** reviewers against the
settled, committed tree — a scientific due-diligence reviewer and a forensic code auditor,
neither given the ledger, the narrative, or any of the reasoning that produced the state
they were reading. That step was the whole point of the exercise, and it worked: they
converged independently on findings that two days of self-audit had missed, including two
that the self-audit had *itself created*.

This section records what they found and what was done, in the same posture as everything
above. The verdicts were "not audit-washing, but incomplete in the direction that flatters
the model" and "fund with conditions (~65%)".

### The two findings the audit created while auditing

**1. A constant was fitted to two benchmarks, and those two benchmarks were then scored as
evidence for it.** Wave H re-derived the methional hydrolysate observability factor
(0.0045 → 0.05623) against `spi_hvp_xylose_120C_PMC9905368` and
`wheat_gluten_hvp_xylose_120C_PMC9905368`, and the resulting agreement on those same two
rows became the *only two hits* in the headline "2 of 11 literature rows inside the CI".
A constant solved from a measurement reproduces that measurement. The honest count for that
headline was **0 of 11**. The circularity was introduced by the wave instruction itself —
*re-derive, then regenerate* — and it survived every existing check, because no check knew
what a fit target was.

**2. The two benchmarks it was fitted to describe an experiment that was never run.** Their
DOI is live, resolves correctly, and its metadata matches — which is why the 225-anchor
citation sweep passed them. But the paper (`10.1007/s10068-022-01194-w`) reacts protein
hydrolysates with **glucose and fructose at pH 7.5 for 90 minutes** and reports only
**relative GC-MS peak areas**. The files claim **xylose and cysteine at pH 6.0 for 30
minutes** with **absolute ppb values** and per-analyte uncertainties, for two compounds —
2-furfurylthiol and 2-methyl-3-furanthiol — that **the paper never mentions at all**. A
paper reporting relative peak areas cannot be the source of an absolute concentration by
any route: there is no conversion, no internal standard, no calibration curve. Those six
numbers have no possible provenance.

This is a different failure mode from the fabrications found in round 1, and the difference
is the lesson: **the citation gate checks DOI *identity*, and no offline gate can check DOI
*content*.** Every anchor in this repository that has not been read against its source is
exposed to exactly this, and the population of such anchors is large.

Both files are now quarantined with the full evidence
([data/benchmarks/quarantined/README.md](data/benchmarks/quarantined/README.md)), the
methional constant is **reverted**, and the re-derivation record is marked **RETRACTED**
with its generator refusing to run. The separate sulfur-barrier refit was checked and
**stands**: it fits against one real, verified benchmark (Hofmann & Schieberle 1998), and
the quarantined files were already on its forbidden list.

### What else they found

- **Back-solved constants labelled as literature anchors.** Three benchmarks scored
  Pearson R = 1.000 and max ratio 1.000 because the observable factors under test had been
  solved from their own measured values — one carried at **17 significant figures**. The
  arithmetic is now written out at the constant: the pea reference "yields" are that
  benchmark's measured ppb divided by a single common scale of 1269; the soy factors are
  the soy benchmark's own values over 838.8. They are relabelled `fitted_to_benchmark`, and
  the summary reports them as **fit-recovery (not predictive)** in a column of their own.
- **Hold-out points that are not measurements.** Four of the eight are constructed: two are
  geometric midpoints of 10–12× bands (carried to 17 significant figures, next to an
  invented `measurement_date` of "publication-year-01-01"), and two are the source's
  odour-activity value multiplied by *this repository's own* unsourced hexanal threshold —
  so they move if we correct our own constant. Every hold-out row now carries its
  `value_provenance`, and the report leads with the split.
- **The flagship mechanism cited the wrong paper.** `10.1021/jf60200a038`, given for van den
  Ouweland & Peer 1975 wherever the MFT route is described, resolves to a gossypol/rat
  feeding study. Corrected to `10.1021/jf60199a045`, and the citation gate was widened to
  scan `src/**`, `docs/**`, `README.md` and `AUDIT.md` — which immediately surfaced a
  fabricated Elsevier DOI in the literature review that had been sitting outside every
  previous sweep.
- **Reported numbers that flattered.** The README's "0.15 → 0.26 dex" understated a
  degradation the artifact recorded as larger; "76 peer-reviewed papers evaluated" is not
  supported by the file it cites (41 entries, 37 scored); the `no_verifiable_source`
  population was under-counted; the ±110× uncertainty band was clamp-limited to ±100×, with
  10.7% of Monte-Carlo draws pinned at the clamp. All corrected, all measured rather than
  estimated.
- **Chemistry and UX defects**: MFT gated on a hydrogen pool only pyrazine chemistry
  supplies; a lipid Schiff-base SMARTS matching 83% of the core network against a comment
  claiming exclusions it never implemented; a demoted legacy shortcut still firing for
  pentoses; broken Arrhenius pairs (an Amadori half-life of 0.13 ms; a lipid homolysis
  kinetically dead on cooking timescales); the off-note lane advertised in the quickstart
  when no fatty acid is registered as a precursor; an experiment card prescribing an SPI/PPI
  matrix protocol for a free-precursor aqueous system; and `--report` hard-crashing on the
  documented conda path for want of `dvipng`.
- **A fabricated showcase.** The generator behind the README's sample table and figures
  built its results from hard-coded numbers and never touched the model. It announced
  *"DECISION CONFIDENCE: HIGH_CONFIDENCE (82/100)"* and *"SCIENTIFIC ENVELOPE: TRUSTED"*
  against a README stating there is no "high" tier, and attributed the flagship thiols to a
  matrix anchor containing no sulfur chemistry. Because the numbers were invented, the
  tracked figures dated from May and had survived every recalibration since — the one
  artifact guaranteed never to show a regression. It now runs the pipeline.
- **A disabled cap kept disabled partly to preserve a headline.** Two docstrings claimed the
  lipid aldehyde-conversion saturation cap was enabled; it ships `null`. The rationale in
  `data/lit/lipid_oxidation_calibration.json` gives a defensible diagnosis (the hold-out
  failures sit in the same progress regime as trace in-panel points, so it is a calibration
  gap rather than a kinetic-shape problem) — and then ends *"The cap is therefore SHIPPED
  DISABLED (null) … so the headline stays 37/48."* That clause is now quoted in the README
  rather than paraphrased away. "So the headline stays" is not a scientific argument, and it
  should not have been part of one.

### What we found while fixing it

Remediation surfaced two defects neither reviewer saw, both in machinery that had been
trusted because it ran without error:

- **70% of the Monte-Carlo barrier channel was inert.** Ten of the fourteen barrier-family
  priors resolved to keys the engine never emits, so perturbing them moved nothing. Every
  published confidence interval was narrower than the priors it claimed to sample. Note the
  direction: this made coverage numbers **too low**, and fixing it *raises* them — which is
  why it is reported as an interval-width change and never as an accuracy improvement.
  Median CI width went 0.99 → 1.15 dex; no prediction moved.
- **Regenerating the hold-out bundles silently deleted their own provenance.** The
  materializer emitted a fixed set of fields, so a regeneration dropped the typed identifier
  of a DOI-less thesis, wrote `source_doi: ""` in its place — which would have reclassified a
  real external hold-out as internal/synthetic — and erased a hand-written citation audit
  note that contradicted a mis-applied DOI repair. The stale artifact had been masking a
  live data defect.

### The new gates

Round 1 added gates for the failure modes it found. Round 2 found that the audit's own
process could manufacture a new one, so
[scripts/ci/fit_target_gate.py](scripts/ci/fit_target_gate.py) now blocks the class
directly. It refuses a build in which a benchmark is both a fit target and a scored row
without disclosure, requires every fit record to state its **leverage** (free parameters per
fitted row), and rejects shipped calibration constants carrying more than six significant
figures — because a 17-digit constant is the residue of a solve, and the precision is the
tell.

The leverage rule is deliberately two-sided. A fit with enough freedom to reproduce its
target row by row is excluded from the coverage count entirely; a single global scale fitted
across two dozen rows is *not*, because excluding it would delete genuine failures rather
than expose them. Applied honestly, this has an uncomfortable consequence worth stating
plainly: **no literature benchmark in the panel is fully out-of-sample** with respect to the
projection scale constant, and the only genuinely out-of-sample evidence in the repository
is the hold-out set.

### The one result that got better

Not everything in Round 2 moved downward, and the exception is instructive. The flagship
compound was gated on a hydrogen pool whose only source, in a ribose/cysteine system, was
pyrazine chemistry — so MFT went to zero whenever the pyrazine lane was off. Sourcing that
reducing equivalent from the thiol redox couple instead (2 cysteine → cystine + 2[H], exactly
atom-balanced) took `cys_ribose_140C_Hofmann1998` from **5.58× under to 1.45× under, and FFT
from 3.26× under to 1.10× over — with no barrier changed.**

*(Round 3 addendum: the route those numbers were measured on has since been retired on
isotope evidence, and the benchmark now reads 2.25× under / 1.22× over. The structural point
above is unaffected — the reducing-equivalent source was and is a real defect fixed — but the
absolute numbers in this paragraph are superseded. See Round 3.)*

That is worth more than the agreement it produced, because of what it invalidates. Wave H had
concluded, from a careful refit, that *no barrier value in any defensible range* could
reproduce those yields and that the residual must therefore be a volatile-budget allocation
deficit. That conclusion was measured against a network in which the route was structurally
starved. A saturating fit profile is exactly what a throttled route looks like from inside a
refit — and the refit could not see the throttle. `thiol_addition_norfuraneol` now ships at a
value fitted against a network that no longer exists, and re-running that refit is an open
item deliberately **not** closed in this wave: a refit stacked on a chemistry change, in one
pass, makes it impossible to say afterwards which of the two produced the agreement.

### What Round 2 cost, in numbers

Every headline in the table above moved in the same direction, and each fall is a fabricated
or circular support being removed rather than a regression:

| | Before Round 2 | After |
| --- | --- | --- |
| Panel size | 16 benchmarks | **14** |
| Literature rows inside 90% CI | 2/11 | **1/3** (fitted rows removed from both sides) |
| Benchmarks without blocking gaps | 7/16, undifferentiated | **0/6 predictive** |
| Genuine-extrapolation hold-out | 2/5 | **0/5** at the pre-widening prior |
| Pentose ≫ hexose ordering | 15.8× | **8.98×**, of which ~1.29× is structural |
| `no_verifiable_source` references | 43 | **84** (62 still feeding runtime numbers) — revised to **102 / 80 / 62 runtime** later the same day when `data/qm/` was tracked for the first time, then to **120 / 98 / 80 runtime** by Wave T3's labelling pass (2026-08-27); see Open items |
| Papers in the literature review | "76 peer-reviewed" | **41 entries, 37 scored** — the claim is withdrawn |
| Test suite | 1178 passed, 1 environmental failure | **1245 passed, 0 failed** |

The last row is the only one that improved, and it improved because `--report` stopped
crashing on the documented install path — a bug that had been carried as "environmental" for
the entire audit while it silently meant the published figures could not be regenerated at
all.

## Round 3: content verification + route correction (2026-08-27)

Round 2 ended on a sentence that turned out to be a prediction: *"the citation gate checks
DOI **identity**, and no offline gate can check DOI **content**. Every anchor in this
repository that has not been read against its source is exposed to exactly this, and the
population of such anchors is large."* Round 3 tested that by reading the papers.

Two independent passes: one that opened every live-panel and hold-out benchmark's actual
source text and compared it line by line against the repository's stored values, and one that
scored the reaction network's **topology** against the isotope-labelling literature — the
class of error no benchmark can catch, because a model can reproduce a concentration through
a mechanism that does not exist.

### Pass 1 — the numbers were wrong about a quarter of the time

**4 fatal + 2 serious findings on ~24 checkable claims: a ~25% content-error rate.** The
access route that made it possible is worth recording, because three of the four fatal
findings were behind it: **Europe PMC `fullTextXML`** returns the complete JATS XML —
including every data table — for MDPI papers that return HTTP 403 to direct requests. Four
sources that had been unverifiable for the whole audit opened on the first try.

| # | Finding | Correction |
| --- | --- | --- |
| 1 | **FATAL** — Pratap-Singh 2021 (Molecules 26:4104, Table 1, PMC8271896): the paper's hexanal is 1138.00 ± 297.30 ppb (pea) and 1621.71 ± 159.69 ppb (soy). The repository had 260 / 380 — **4.38× and 4.27× low**, origin unknown. Not a units error: the same files' 2-pentylfuran matched the paper exactly. | Corrected; ranks re-derived |
| 2 | **FATAL** — same two files: hexanol 80 / 120 ppb was **fabricated**. Table 1 says *n.d.* for both, and the text says pea proteins "contained no alcohol compounds". Both files ranked hexanol at `expected_rank 3`. | Rows removed |
| 3 | **FATAL** — Li 2026 hold-out: two of four points were **the wrong table rows**. Nonanal 29.42 was the *Decanal* row (true 72.66); 2-pentylfuran 221.5 was the *Maltol* row (true **5625.80**, 25× low). Both were labelled `reported_point_value`. | Corrected against PMC12984281 |
| 4 | **FATAL** — Ma 2024 acrylamide 62.62 ppb **is not in the paper**. Acrylamide appears only in Figure 2D, whose 130 °C bar reads ~150 µg/kg. Every plausible derivation from the companion Fu 2023 was ruled out. The benchmark id says "ACSRef3" against an MDPI DOI — the value predates the DOI and was never re-checked. | Re-provenanced as `figure_readoff`, 150 ppb, tolerance widened 15% → 20% |
| 5 | **SERIOUS** — Hernandez/Resconi furfural 1040 ppb was the mean of **2 of 3** products, silently dropping Impossible Burger (64.71). A cherry-pick sitting inside a ±5% tolerance. | Honest 3-product mean 715.22, uncertainty 5% → 79% (the 17× inter-product spread is real) |
| 6 | **SERIOUS → FATAL, escalated 2026-08-27 (Waves S2b/S2c).** Originally read: *"Hofmann & Schieberle 1998 reports **mol %**, not ppb. The 342 / 200 ppb values require a conversion documented nowhere in this repository, on the panel's tightest contract (1.45×)."* The escalation: there is nothing on the far end of the conversion. Wave S2b traced 342 and 200 to `data/benchmarks/maillard_validation_benchmarks.md` §1.3 — an abstract-reconstructed range table committed in the **same commit** as the benchmark — and both are interior points of two invented, overlapping mol % bands. Not a units problem; a fabrication. | Round 3: values unchanged, `content_verification_note` added. **Wave S2c:** both values marked `no_verifiable_source` (unedited), the 1.45× / 0.09 dex contract **retired**, tier demoted `PRIMARY → REFERENCE`, `thiol_addition_pentodiulose` reverted 26.35 → 28.60 and its fit record retracted. §1.3 kept under a fabrication warning |

Nothing about these was detectable offline. Every DOI resolved, every metadata field matched,
and the citation gate passed all six.

### Pass 2 — the flagship mechanism was the wrong mechanism

28 model routes and route-absences were scored against isotope-labelling and
pathway-elucidation studies: 8 confirmed, 4 contradicted, 10 partially supported, 3 untested,
1 correct omission, 3 completeness gaps. Full dossier:
[docs/validation/isotope_topology_evidence.md](docs/validation/isotope_topology_evidence.md).

The contradicted one is the flagship. Round 1 had replaced a fabricated one-step MFT shortcut
with the "accepted" **1-deoxyosone → norfuraneol → +H₂S → MFT** route, cited to van den
Ouweland & Peer 1975 — and Round 2 had even repaired that citation from a wrong DOI to the
right one. The paper is real, the citation is now correct, and **the route is still wrong**:
van den Ouweland & Peer describes a *synthesis* of MFT from norfuraneol, not the in-situ
Maillard pathway.

- **Cerny & Davidek 2003** (`10.1021/jf026123f`) spiked authentic unlabelled norfuraneol into
  a [¹³C₅]ribose/cysteine system. The MFT came out **mainly ¹³C₅-labelled** — from the ribose,
  not from the spike: "4-hydroxy-5-methyl-3(2H)-furanone is **unimportant as an
  intermediate**". They propose 1,4-dideoxypento-2,3-diulose instead.
- **Cerny & Davidek 2004** (`10.1021/jf035265m`) confirms that intermediate positionally with
  [1-¹³C]ribose.
- **Hofmann & Schieberle 1998** (`10.1021/jf9705983`) — the paper this repository's only
  surviving sulfur anchor comes from, and the one its sulfur barrier was fitted against —
  says in its own abstract that norfuraneol/cysteine is among the *less* effective MFT
  precursors. That had been sitting in the repository, uncontradicted, the whole time.

`Thiol_Addition_Norfuraneol` was retired. Two families replace it: `Deoxyosone_Reduction`
(1-deoxy-2,3-pentodiulose + 2[H] → 1,4-dideoxypento-2,3-diulose + H₂O, barrier 28.0
ESTIMATED) and `Thiol_Addition_Pentodiulose` (1,4-dideoxyosone + H₂S → MFT + 2 H₂O, barrier
28.60 ESTIMATED and explicitly **unconstrained**). Both are atom-balanced and RDKit-verified.
The sulfur-incorporation step now needs **no** reducing-equivalent token at all — the
literature topology moves that bookkeeping upstream onto a dehydration where a reductone
donor is chemically ordinary, so the correction *removes* a lumping fiction rather than adding
one. Norfuraneol is kept as a genuine furanone product; it just no longer feeds MFT.

### What Round 3 cost, in numbers

| | Before Round 3 | After | Cause |
| --- | --- | --- | --- |
| Hofmann MFT | 235.32 ppb — **1.45× under** | 151.87 ppb — **2.25× under** | Route correction |
| Hofmann FFT | 219.96 ppb — 1.10× over | 243.72 ppb — **1.22× over** | Route correction |
| Pentose ≫ hexose ordering | 8.98× (ribose 981.3) | **3.39×** (ribose 370.3), of which ~1.13× structural | Route correction |
| Benchmarks without blocking gaps | 0/6 predictive, **3/4** fit-recovery, 4/4 synthetic (7/14) | 0/6 predictive, **1/4** fit-recovery, 4/4 synthetic (**5/14**) | Pratap-Singh content correction |
| Pea / soy ambient slurry | max ratio 1.002× / 1.001× (`pass`) | **4.37× / 4.27× under** (`scale-gap`) | Pratap-Singh content correction |
| Pea / soy ambient slurry, after the Wave O refit | 4.37× / 4.27× under (`scale-gap`) | **1.0113× on both** (`pass-no-ranking`, still `fit_recovery`) | Observability refit to the verified anchors |
| Acrylamide extrusion | 6.42× under | **15.39× under** | Ma 2024 re-provenance (62.62 → 150 ppb) |
| Resconi furfural | 3.14× over | **4.56× over** | Cherry-pick corrected (1040 → 715.22 ppb) |
| Hold-out median fold error | 32.79× | 15.31× → 42.62× after Wave O → **93.68×** after Wave R | Li 2026 wrong-row correction, then the Wave O observability refit, then the Wave R Liu reference correction |
| Hold-out genuine extrapolations | 0/5 pre-widening | **1/5** pre-widening (unmoved by Wave O) | Li 2026 wrong-row correction |
| Hold-out coverage at the shipped sigma | 5/8 | 4/8 after Wave O → **3/8** after Wave R | Wave O observability refit, then the Wave R Liu reference correction |
| Literature-row projection objective | 0.74 dex | 0.96 dex → **0.88 dex** after Wave O | Both corrections; then the Wave O refit (the 0.08 dex gain is entirely two fit-recovery rows and is **not** evidence) |
| Matrix leave-lane-out ln-sigma | 2.86 RMS on 6 residuals | **3.25** on 5 | Pratap-Singh content correction |
| Matrix leave-lane-out ln-sigma, after the Wave O refit | 3.25 on 5 | **3.25 on 5, bit-identical** | The derivation never reads the observability factors — by design |
| Test suite | 1245 passed, 0 failed | **1242 passed, 1 skipped, 2 xfailed, 0 failed** | 18 stale pins re-derived, none relaxed |
| *(Wave P)* Test suite | 1242 passed, 0 failed | **1265 passed, 1 skipped, 2 xfailed, 0 failed** (873.99 s) | 1265+1+2 = 1268 collected = Wave O's 1248 + the 20 new Wave P tests. One real failure was produced and fixed en route: the Trikusuma nonanal recovery, a stale pin from the oleate substrate correction, re-pinned two-sided. Three earlier attempts were abandoned to a machine condition, not the tree — the data volume hit 100% (ENOSPC produced spurious ERRORs that do not reproduce) and free RAM fell to ~60 MB. |
| *(Wave P)* Hofmann MFT | 151.87 ppb — 2.25× under | **242.38 ppb — 1.41× under** | `thiol_addition_pentodiulose` refit 28.60 → 26.35 (fit recovery, not evidence) |
| *(Wave P)* Hofmann FFT | 243.72 ppb — 1.22× over | **217.99 ppb — 1.09× over** | Same refit; FFT was **not** fitted and co-moves through shared upstream flux |
| *(Wave P)* Pentose ≫ hexose ordering | 3.39× (1.13× structural) | **6.15×** (2.31× structural) | The refitted barrier gap, not the mechanism |
| *(Wave P)* Benchmarks without blocking gaps | 0/6 predictive, 1/4 fit-recovery, 4/4 synthetic (5/14) | 0/6 predictive, **0/4** fit-recovery, 4/4 synthetic (**4/14**) | Oleate substrate correction broke the last fit recovery |
| *(Wave P)* Trikusuma heated pea nonanal | 24.00 ppb (`pass`, exact recovery) | **10.56 ppb — 2.2727× under** (`scale-gap`) | Nonanal moved onto the oleate pool; 2.2727 = 1/(22.0/50.0) exactly |
| *(Wave P)* Hold-out Li 2026 nonanal | 272.63× over | **118.31× over** | Oleate substrate correction — nothing fitted |
| *(Wave P)* Hold-out Liu 2023 nonanal | 10.86× over | **4.78× over** | Same |
| *(Wave P)* Hold-out headline | median 42.62×, coverage 4/8, worst 2474× | **all three unchanged** | The median sits between two points that did not move |
| *(Wave R)* Hold-out headline | median 42.62×, coverage 4/8 | **median 93.68×, coverage 3/8** | The Liu bundle's two reference values were replaced with the ones in its source; no prediction moved |
| *(Wave P)* Matrix leave-lane-out ln-sigma | 3.25 on 5 residuals | **3.02** on 5, 90% CI [2.03, 6.30] | The Trikusuma nonanal residual fell 3238.93 → 1425.13 ppb |
| *(Wave P)* Literature-row projection objective | 0.88 dex | **0.89 dex** at the shipped tau; the fit optimum moved 5011.87 → 10000 min (2.51× → **1.26×** away) | Chemistry additions + the sulfur refit |
| *(Wave S1)* Hofmann MFT | 242.38 ppb — 1.4110× under | **283.59 ppb — 1.2060× under** | Additive flux propagator — no constant moved |
| *(Wave S1)* Hofmann FFT | 217.99 ppb — 1.0900× over | **297.28 ppb — 1.4864× over** | Same; **worse**, and not clawed back |
| *(Wave S1)* Hofmann contract | max ratio 1.4110 / MALE 0.0935 (fails on **one** criterion) | **1.4864 / 0.1267 dex** (fails on **both**) | Same — the untouched contract is 1.45× / 0.09 dex |
| *(Wave S1)* Pentose ≫ hexose ordering | 6.15× (2.31× structural) | **7.78×** (**3.14×** structural) | The pentose limb has more parallel routes; the propagator stopped discarding them |
| *(Wave S1)* Cerny 2008 MFT (reference-only) | 3.195× over | **2.787× over** | Additive flux propagator |
| *(Wave S1)* Bolton 1994 MFT | 763.588× under | **748.022× under** | Same — a 2% move on a 750× miss changes nothing |
| *(Wave S1)* Resconi furfural | 4.657× over | **4.402× over** | Same — furfural loses budget share to the sulfur channels |
| *(Wave S1)* Internal snapshots, Hexanal | pea 0.1720 / soy 0.1783 ppb | **pea 0.7425 (×4.317) / soy 1.7006 (×9.540)** | The matrix calibration registry became reachable on the augmented lane; these are the Wave O refitted factors finally applying |
| *(Wave S1)* Hold-out headline | median 93.68×, coverage 3/8, worst 2474× | **all eight points bit-identical** | Every hold-out bundle runs the `matrix_only` path, which neither fix touches |
| *(Wave S1)* Matrix leave-lane-out ln-sigma | 3.02 on 5 residuals | **bit-identical** | The derivation reads neither the propagator nor the observability factors |
| *(Wave S1)* MC external-literature CI width | 0.8495 dex | **0.9463 dex**, coverage still 1/3 | Additive channels sample more barriers, so intervals widened **without** buying coverage |
| *(Wave S1)* Benchmarks without blocking gaps | 0/6 predictive, 0/4 fit-recovery, 4/4 synthetic (4/14) | **unchanged** | The four synthetic snapshots were refreshed, as they are after every intentional model change |
| *(Wave S1)* Mechanistic-priority benchmarks | 2 | **0** | A side effect of the registry repair, and a **loss of signal**: `_matrix_closure_action` has no branch for `process_state_mismatch`, so a factor reached only through a state fallback now scores as an acceptable class-level transfer. Open item. |
| *(Wave S1b)* **Directional panel, pH claims** | **2/7** | **4/7** | The pyrazine ionisation branch and the enolisation route-selection term were both disconnected; both are now wired |
| *(Wave S1b)* Directional panel, moisture claims | 0/3 | **0/3** | The aw correction now reaches the chemistry, but HMF — the observable in 2 of the 3 rows — shares its producing family with furfural, so their shares are pinned against each other |
| *(Wave S1b)* Directional panel, headline | 19/29 (66%) | **20/29 (69%)** | +2 on pH (PH-04, PH-06), −1 on sugar identity (SUG-12) |
| *(Wave S1b)* Directional panel, non-pH/aw | 17/19 (89%) | **16/19 (84%)** | SUG-12: glucose's HMF route gets the acid boost, fructose's direct dehydration is excluded from that set |
| *(Wave S1b)* Hofmann MFT | 283.59 ppb — 1.206× under | **154.85 ppb — 2.209× under** | The 1,2-enolisation acid boost at pH 5.0 moves budget share from the MFT arm to the FFT/furfural arm. **Worse; not clawed back** |
| *(Wave S1b)* Hofmann FFT | 297.28 ppb — 1.486× over | **267.50 ppb — 1.338× over** | Same mechanism, opposite sign |
| *(Wave S1b)* Hofmann contract | 1.4864 / 0.1267 dex (fails on both) | **2.2086 / 0.2352 dex** (fails on both, by more) | Untouched contract is 1.45× / 0.09 dex |
| *(Wave S1b)* Bolton 1994 MFT | 748.02× under | **6730.85× under** | Almost entirely the aw term: at aw 0.98 the thiamine→MFT lane's water-shedding aromatisation takes a ×0.32 penalty its upstream steps do not |
| *(Wave S1b)* Cerny 2008 MFT (reference-only) | 2.787× over | **23.406× over** | Same |
| *(Wave S1b)* Resconi furfural | 4.402× over | **5.457× over** | Furfural gains share from the acid-favoured 1,2-enolisation route |
| *(Wave S1b)* Pentose ≫ hexose ordering | 7.78× (3.14× structural) | **18.27×** (**4.27×** structural) | **Not better discrimination.** Both numbers fell (824.7 → 374.0, 106.0 → 20.5); the ratio rose because the hexose limb's three-water shortcut is penalised harder at aw 0.98 |
| *(Wave S2c)* Pentose ≫ hexose ordering | 18.27× (4.27× structural) | **8.26×** (**4.27×** structural) | The Wave P barrier refit was reverted as circular. Only the pentose limb runs it: ribose 374.0 → 169.1 ppb, **glucose unchanged at 20.5**. Smaller claim, better evidence — the mechanism-carried share went 23% → 52% |
| *(Wave S2c)* Hofmann's status | the panel's tightest literature contract (1.45× / 0.09 dex), `tier: PRIMARY` | **not literature at all**; contract **RETIRED**, `tier: REFERENCE`, both values `no_verifiable_source` | 342 = 0.0300 mol % × 0.010 M × 114.17 and 200 = the geometric mean of an invented `~0.01–0.03` band. The contract was ~1.7× tighter than the 2.5× spread of the band it came from |
| *(Wave S2c)* `thiol_addition_pentodiulose` | 26.35 kcal/mol, **fitted** | **28.60 kcal/mol, estimated** | Its sole fit target is now known non-evidence. Record retracted; same treatment Wave I gave the Methional re-derivation |
| *(Wave S2c)* Hofmann MFT | 154.85 ppb — 2.2086× under | **78.09 ppb — 4.3797× under** | The barrier revert. **Much worse, and published worse** — but the yardstick is fabricated, so it is not evidence about the chemistry |
| *(Wave S2c)* Hofmann FFT | 267.50 ppb — 1.3375× over | **293.67 ppb — 1.4684× over** | Same; the two lanes share their upstream trunk |
| *(Wave S2c)* Hofmann contract | 2.2086 / 0.2352 dex against 1.45× / 0.09 dex | **4.3797 / 0.4041 dex** against the inherited global default 1.5× / 0.10 dex | Retired, **not widened**. The fallback default is marginally looser and is failed by more; that inheritance is stated in the file rather than left to be found |
| *(Wave S2c)* Cerny 2008 MFT (reference-only) | 23.406× over | **25.741× over** | Same barrier; the xylose lane runs the same pentodiulose step |
| *(Wave S2c)* Strict-ready | 0/14 | **0/14** | Retiring a *failing* contract removes a failure, not a pass |
| *(Wave S2c)* MC honest literature coverage | 0/3 | **0/3** | The benchmark stays a declared fit target of two other live records, so it stays out of the numerator **and** denominator — where a value with no verifiable source belongs |
| *(Wave S2c)* Hold-out headline | median 93.68×, coverage 3/8 | **all eight points bit-identical, third wave running** | `matrix_only` never reaches `predict_from_steps`, so no barrier is on its path |
| *(Wave S2c)* MFT lane rate-limiting step | `Amadori_Rearrangement` (shared with the FFT lane) | **`Thiol_Addition_Pentodiulose`** (29.5354 effective vs Amadori's 29.0603) | A real structural movement found by a re-pinned test. It does not reopen the propagator rule — channel ids are full step-sets — but it retires the "both lanes bottleneck on the same trunk" argument |
| *(Wave S2c)* Internal snapshots | pass 4/4 after refresh | **pass 4/4 after refresh** | MFT ×0.416079 and the MFT-derived disulfide ×0.420049; **every** other propagated row ×1.010812 — one shared factor across three unrelated compounds is the fixed budget re-normalising |
| *(Wave S1b)* Internal snapshot ranking | MFT > FFT, DMP > disulfide | **FFT > MFT, disulfide > DMP** | A real ranking movement, not a refresh artifact |
| *(Wave S1b)* Benchmark `status_counts` | scale-gap 8 / pass-no-ranking 2 / pass 4 | **unchanged** | Four rows degraded within the scale-gap bucket; the four synthetic snapshots recover after the documented refresh |

**The Hofmann row is the one to read twice.** The 1.45× was bought with
`thiol_addition_norfuraneol` = 26.85 kcal/mol, a barrier fitted *through* the route the
isotope evidence contradicts. Its replacement ships at the un-fitted sulfur-addition class
value. Nothing was tuned to recover the lost agreement, and nothing should be: **a worse
number obtained through a route the literature supports is worth more than a better number
obtained through one it refutes.**

**And the fit-recovery row is the sharpest result in three rounds.** Round 2 relabelled the
Pratap-Singh lanes `fitted_to_benchmark` and reported their 1.002× / 1.001× agreement as
"fit-recovery (not predictive)" rather than as evidence. Round 3 shows that even that was
generous: the constants were back-solved from values the paper does not contain, so with the
paper's real numbers they miss by exactly the size of the transcription error. **A
fit-recovery that no longer recovers is proof that the recovery was never evidence.** The
observability factors were deliberately *not* refitted — doing so in the same pass as a
chemistry change would make the resulting agreement unattributable — so
`src/matrix_calibration_registry.py` now carries the arithmetic and the consequence at the
constant.

### Wave O — the refit, and the price it charged (2026-08-27, owner-approved)

The deferred refit was then done, on its own, against the verified anchors.

**What was fitted.** One parameter: a shared multiplicative scale of **4.317249×** on the
ambient-slurry *hexanal* observability factors of both lanes (pea `1.0 → 4.31725`, soy
`2.2097561 → 9.54007`; the soy heated-matrix factor is *defined* as the soy ambient baseline
× the Shu 2024 attenuation, so it was propagated `0.649668 → 2.80478` rather than left
composed with a stale baseline). Nothing else moved — no barrier, no projection constant, no
marker yield, and not the 2-pentylfuran factors, whose anchors (638 / 2492 ppb) Wave K
verified verbatim. Generator and record:
`scripts/generators/refit_matrix_observability_pratap_singh.py` →
`results/validation/matrix_observability_refit_pratap_singh.{json,md}`.

**One parameter, not two, deliberately.** Two free factors against two rows is an exactly
determined system: its zero residual is arithmetic and measures nothing. One shared scale
leaves the data a degree of freedom to disagree with, and the residual it reports is
informative: **1.0113× on both rows**. The pea lane wanted 4.36606× and the soy lane
4.26899×, so the two corrected anchors are mutually consistent to 1.1%. The transcription
error was a common *absolute-scale* error; the pea-vs-soy release ratio the registry encodes
(2.2098) survived it untouched. Nothing saturated and no bound was hit — the optimum sits
more than a decade inside the search range on both sides.

**The status did not change, and that is the point.** Both benchmarks were already
`fitted_to_benchmark`; they still are, still sit in the `fit_recovery` bucket, and are still
excluded from the honest coverage numerator *and* denominator. Their `overall_status` moved
`scale-gap` → `pass-no-ranking`, which is **not** counted as a pass: 0/6 predictive and
5/14 aggregate are unchanged. A refit changes the size of a recovery, never its evidential
status.

**The price, which is the number that matters.** The external hold-out — the only lane a
refit cannot flatter — got **worse**:

| Hold-out point | Before | After |
| --- | --- | --- |
| Bi 2020 raw pea, hexanal | 5.37× | **1.24×** |
| Liu 2023 PPI, hexanal | 4.52× | **19.50×** |
| Li 2026 HME, hexanal | 21.58× | **93.15×** |
| the other five points | unchanged | unchanged |
| **median fold error** | 15.31× | **42.62×** |
| coverage at the shipped sigma | 5/8 | **4/8** |
| worst miss / pre-widening 1/5 | 2474× / 1 of 5 | unchanged |

The cause is a contradiction inside the literature, not inside the fit. The pea ambient lane
carries two external measurements of nominally the same system that disagree by **24×** —
Bi 2020 at 1260 ppb and Liu 2023's band midpoint at 51.96 ppb — and the erroneous 260 ppb the
old constants reproduced sat almost exactly at their geometric mean (√(1260 × 51.96) = 255.9).
**Being wrong in the middle of a contradiction scores better than being right at one end of
it.** The verified anchor (1138 ppb) agrees with Bi to 1.11× and sits 6.3× above the top of
Liu's reported band; which is representative of commercial PPI is a question this repository
cannot settle. What the refit bought is that the shipped constants are now anchored to
numbers that exist in a paper.

> **Superseded 2026-08-27 (Wave R): the contradiction was not in the literature either.**
> The "Liu 2023 band" in the paragraph above appears in no source. The thesis was retrieved
> and read in full and its Table 2.7 reports hexanal at 2445–52454 µg/L, so the verified
> 1138 ppb anchor sits *just under Liu's lowest lot* rather than 6.3× above her band, and Bi
> and Liu do not disagree by 24× on any comparable basis. The Wave O refit's direction is
> **vindicated by the corrected target**: the same Liu hexanal point reads 19.50× → **11.17×**
> once scored against the real number. What still costs the headline is the *nonanal* row,
> not a clash between two believable papers — see "Wave R" below.

**What the refit could not touch.** The 1-hexanol factors. The paper reports n.d. in both
matrices, so the shipped soy value `0.143 / 0.063` is a ratio of two numbers (120 and 80 ppb)
that appear nowhere in it, and there is no anchor to refit it against. It is left in place,
flagged in the registry and in the fit record, and it is the constant behind the hold-out's
1117× 1-hexanol miss on Li 2026. Retiring it is a separate science decision.

**And the sigma did not move at all.** Re-deriving the leave-lane-out matrix ln-sigma after
the refit produced a **bit-identical** artifact (RMS 3.2520 on 5 residuals, 90% CI
[2.186, 6.796], shipped 2.86 still inside). That is structural: the derivation multiplies the
oxidation load by the base marker yield and never reads an observability factor, because the
uncalibrated tier exists to describe a lane that has no such calibration. **No refit of these
constants can ever be used to justify narrowing that prior** — now stated in the generator's
own header.

### Wave P — evidence-grounded chemistry additions (2026-08-27, owner-approved)

Six changes, each anchored to a citation retrieved and CrossRef-verified during the wave,
none tuned to reproduce a prior output. Full per-item record: `tasks/audit_remediation.md`,
Wave P.

1. **`thiol_addition_pentodiulose` refitted** against `cys_ribose_140C_Hofmann1998` alone,
   28.60 → **26.35** kcal/mol. One free parameter; range [23.30, 29.65] bounded by values
   already in `FAST_BARRIERS`; decision rules copied from the Wave H script; run **last** in
   the wave so the fit sees the network that ships. Objective 0.2192 → **0.0935 dex**. The
   profile minimum sits at the range **floor** (`hit_a_bound = true`), so the adopted value is
   the conservative edge and the residual is not removable by this barrier. The fit target's
   own unverified mol%→ppb conversion now travels **verbatim** inside the constant's rationale,
   which localises that risk in one named number. Record:
   `results/validation/sulfur_barrier_refit_pentodiulose.{json,md}`, superseding the stale
   `sulfur_barrier_refit_hofmann.{json,md}`.
2. **The C2 + C3 recombination lane to MFT** — Hofmann & Schieberle 1998's own
   highest-yielding system (1.4 mol %, 6 min, 180 °C, dry), which the network could not
   express because mercapto-2-propanone was not a species. Three balanced steps, mechanism
   from Cerny 2015 (`10.1016/b978-1-78242-103-0.00009-6`, full text) rendering their scheme.
   **Its measured contribution to the shipped MFT number is exactly zero** — see the finding
   below.
3. **Norfuraneol's real sulfur fate.** Wave N retired it from the MFT lane and left it with no
   consumers at all; Cerny & Davidek 2003 (`10.1021/jf026123f`) assign it 2-mercapto-3-
   pentanone, and Whitfield & Mottram 1999 (`10.1021/jf980980v`) supply the 2,3-pentanedione
   intermediate they measure in exactly that system.
4. **Nonanal is cleaved from the oleate pool**, not the linoleate pool. `oleic_acid_pct` was
   dead code. Two hold-out points improved with nothing fitted; the last fit-recovery pass in
   the panel broke.
5. **Fructose reaches HMF by its own ring-retained dehydration**, not through the glucose
   3-deoxyosone that Perez Locas & Yaylayan 2008 (`10.1021/jf8010245`) label-exclude for it.
   Zero live benchmark contains fructose, so this moves no scored row — it is a topology
   repair, and it is reported as one.
6. **The furaneol `[HH]` pool gate is gone.** Red-team finding H4's second half: predicted
   furaneol from glucose was contingent on pyrazine chemistry. The token was not re-sourced,
   it was removed, and the amino acid was written into the step as the reducing partner.
   Measured: disabling aminoketone condensation now leaves the DMHF step standing (it used to
   take it to zero).
   **Citation corrected 2026-08-29 (Wave Q1).** This item used to read "the accepted mechanism
   names the reductant and it is the amino acid (Blank & Fay 1996 …)". Blank & Fay 1996 do **not**
   name it: they leave it open — *"The reduction may occur either by a dismutation or by a
   reaction with further enoloxo compounds"* (`10.1021/jf950439o` p. 534, after Schieberle 1992),
   under a figure legend their own paper calls *"Hypothetical"*. In that paper the amino acid is
   the **carbon** donor (the Strecker aldehyde), not a hydride donor. The coupling the model
   ships rests on **Kerler et al. 2010** (`10.1002/9781444317770.ch3`) alone, and on a weaker
   statement than "names the reductant" — Furaneol from acetylformoin *"was significantly
   enhanced in the presence of reductones such as ascorbic acid or methylene reductinic acid as
   well as the Strecker-active amino acid proline"*. It is therefore a **declared modelling
   choice**, not a mechanism a cited paper asserts. Resolving it properly needs Hofmann &
   Schieberle 2001, *"Acetylformoin — an important progenitor"*, Flavour 2000 Proc. 311–322,
   which this repo does not hold.

**The finding this wave produced, and it is about the engine rather than the chemistry.**
Adding a second, literature-evidenced channel to the flagship compound changed its prediction
by **exactly nothing**. Measured on Hofmann1998 at the shipped constants: the pentodiulose lane
alone gives 242.38 ppb, the new C2 + C3 lane alone gives 71.02 ppb, and **both together give
242.38 ppb**. `src/recommend.py` relaxes to the *lowest-span path per product*; it does not sum
flux over parallel channels. Had the two added, MFT would read 313.39 ppb (1.09× under instead
of 1.41× under). Every claim this repository makes about a compound therefore rests on one
route, and enriching the network cannot improve a number unless the new route is *faster* than
the old one. Reported, not fixed: making the propagator additive is a model-wide recalibration
event and belongs to its own wave. Pinned in
`tests/scientific/test_wave_p_chemistry_2026_08.py`.
**→ CLOSED by Wave S1 (below). The propagator is now additive, MFT reads 283.59 ppb, and the
313.39 ppb projected in the paragraph above was never obtainable — the two lanes share their
rate-limiting step and the volatile budget is fixed. That paragraph is kept as the record of
the finding; its arithmetic is superseded.**

**And a naming hazard, found while placing the new families.** `src/conditions.py` selects the
pH-ionisation and Labuza water-activity corrections by **substring match on the family name**.
A family called `Furanone_Strecker_Reduction` would have silently collected both (~480×
suppression at pH 5.5 / a_w 0.95) purely because of how it was spelled. The Wave P families are
named to receive the same treatment as the sibling steps they extend, and a test now pins that
none of them picks up a correction. The underlying coupling is untouched and remains open.

### Wave R — the hold-out was being graded against numbers from nowhere (2026-08-27, owner-adjudicated)

Wave O closed with an open `[P]`: *"Bi 2020 (1260 ppb) vs Liu 2023 (15–180 ppb band) disagree
by 24× on nominally the same system and no observability factor satisfies both."* That was
treated for three waves as the sharpest unresolved conflict in the pea lane, quoted in this
document, in the README and in the validation contract as evidence that the literature itself
was inconsistent.

**There was no conflict. One side of it was not in the literature.**

The source was located and retrieved in full: Yaozheng Liu, *Flavor Chemistry of Pea Proteins*,
MS thesis, North Carolina State University, **2021** — published as Liu, Cadwallader & Drake
(2023), *Food Chemistry* 406:134998, `10.1016/j.foodchem.2022.134998` (CrossRef-verified). Its
Table 2.7 quantifies nine commercial pea proteins, rehydrated to 10 % solids, by
HS-SPME-GC-MS/MS against five-point external standard curves.

| Compound | This repository carried | The thesis reports (Table 2.7) | Error |
| --- | --- | --- | --- |
| hexanal | band 15–180 ppb, mid 51.96 | **2445–52454 µg/L** (4318, 3360, 2445, 6052, 6383, 11203, 12181, 52454, 2533) | 50–300× **low** |
| nonanal | band 5–50 ppb, mid 15.81 | **0.188–3.42 µg/L** | 6–266× **high** |
| (E,E)-2,4-heptadienal | band 0.5–8 ppb, OAV 7–114 | **the compound does not appear in the thesis at all** | unsourceable |
| 3-isobutyl-2-methoxypyrazine | max 0.08 ppb | **the compound does not appear in the thesis at all**; its methoxypyrazine is IPMP at 6.126–57.0 µg/L | wrong compound, and 713× low for the nearest real one |
| hexanal OAV (`computational_priors`) | 28, with 1-octen-3-ol at 14 | **543–11656** (Table 2.8), 1-octen-3-ol **2–47** | 19–416× low |

Not one of those numbers is a rounding, a unit slip or a transcription off an adjacent row —
the failure mode Round 3 found in the Li 2026 and Pratap-Singh bundles. They correspond to
nothing in the document they cite. The two compounds that are simply *absent* were **retired,
not repaired**: `no_verifiable_source`, removed from any scored target, with the nearest real
rows recorded in the correction note and explicitly *not* substituted.

**The price, again on the only lane a correction cannot flatter.**

| Hold-out point | Measured, before → after | Fold error, before → after |
| --- | --- | --- |
| Liu 2023 PPI, hexanal | 51.96 → **11320 ppb** | 19.50× → **11.17×** (better) |
| Liu 2023 PPI, nonanal | 15.81 → **0.8018 ppb** | 4.78× → **94.22×** (worse) |
| the other six points | unchanged | unchanged |
| **median fold error** | | 42.62× → **93.68×** |
| coverage at the shipped sigma | | 4/8 → **3/8** |
| worst miss / pre-widening 1/5 | | 2474× / 1 of 5, unchanged |

**No prediction moved.** Nothing was fitted, no constant was touched, and the hold-out stayed
outside every fit — the three gates were re-run green afterwards. The headline got worse
because the grading key was replaced with the real one.

Two things worth stating plainly. First, the **Wave O refit is vindicated by the corrected
target**: the Liu hexanal point improved from 19.50× to 11.17×, and the verified 1138 ppb
Pratap-Singh anchor that drove the refit now sits *just under Liu's lowest lot* instead of
6.3× above her band. The refit was right and was being graded by a broken key. Second, the
**nonanal row is now the sharpest lipid-lane over-prediction the repository has against a
directly-quantified reference**: 75.5 ppb predicted against a band whose top is 3.42 ppb.
Wave P's oleate-substrate correction is the partial mitigation already landed — it took this
same point from 214× to 94× — and it was not enough.

And the caveat travels with the corrected values: Liu's standard curves were built **in
deionized water, not in the protein matrix**, so protein binding of the analyte is
uncorrected and 0.188–3.42 µg/L is a *lower* bound on total nonanal. The over-prediction is,
if anything, understated by this table.

### Wave S1 — the additive flux propagator, and the registry nobody could reach (2026-08-27, owner-approved)

Two **structural** fixes to `src/recommend.py`. No barrier, no observability factor, no
projection constant was touched. Both were findings carried forward from earlier waves as
`[P]`, both were reported-not-fixed at the time, and both are now closed.

**Fix 1 — parallel channels can add.** Wave P's own headline finding: `predict_from_steps`
relaxed to the *lowest-span route per product* and kept only that route's flux, so adding a
real, literature-evidenced second route to the flagship compound contributed **exactly zero**
(pentodiulose lane alone 242.38 ppb, C2+C3 lane alone 71.02 ppb, both together 242.38 ppb).
Every compound-level claim in the repository rested on one route, and *"the model has N routes
to X"* was not a statement about the prediction of X. Each product's flux proxy is now the
**sum** over kinetically distinct routes, deduplicated on the route's full ordered step-set.

> **The obvious guard was implemented, measured, and rejected — and that is the load-bearing
> result.** The natural dedupe rule is *"two routes sharing a rate-limiting step are one
> channel; take the max."* Measured: **both** MFT routes here have the same highest-barrier
> step, the trunk `Amadori_Rearrangement` at 29.06 kcal/mol, which sits on the shared
> cysteine/ribose stem that essentially every route in the network passes through. Under that
> rule MFT keeps its old 242.38 ppb *exactly* and the whole live panel moves 3 rows by less
> than 1.15× — winner-takes-all in all but name. It is also wrong on the physics: for a slow
> trunk feeding two fast branches to P and one to Q, the trunk fixes the total flux and the
> branches **partition** it by conductance, so P's share is 2/3, not 1/2. Because this
> propagator's per-route weight is `pool · exp(−span/RT)` with
> `exp(span/RT) = Σᵢ exp(barrierᵢ/RT)`, a dominant trunk collapses every route's weight onto
> the same trunk value — and it is precisely **summing** them that reproduces the 2/3.

This also **corrects Wave P's arithmetic**, which projected 242.38 + 71.02 = 313.39 ppb "if
the two channels are genuinely independent". They are not independent, and the volatile budget
is fixed, so the single-lane numbers were never addable. The shipped answer is 283.59 ppb.

**Mass honesty, verified rather than asserted.** The allocation layer normalises activities to
mole fractions of `total_volatile_budget_molar`. Measured on three systems (Hofmann, Resconi,
the soy snapshot), running the same converged state through the additive proxy and through the
retained winner-takes-all proxy: **allocated molar / budget = 1.000000000000 in every case,
before and after.** The summed *ppb* moves by ≤0.6 %, and only because ppb is
molecular-weight-weighted and the allocation shifted between species of different mass. Adding
channels moves the split; it cannot mint mass.

**Fix 2 — the matrix calibration registry was unreachable.** Wave O's finding (f): on the
`matrix_precursor_augmented` lane, species injected from the lipid-oxidation path arrive
labelled by *canonical SMILES* (`CCCCCC=O`), so the name-keyed
`describe_matrix_calibration` matched nothing and the observability factor silently applied as
**1.0**. Wave O proved it by changing the soy hexanal factor 4.32× and finding the snapshot
bit-identical. The lookup now resolves a SMILES label to its compound-database name at the
boundary, and the same perturbation moves the prediction by exactly 4.32× — pinned as a test.
The internal snapshots can, for the first time, detect drift in
`src/matrix_calibration_registry.py`.

**What each fix moved, separately attributed.**

| | Fix 2 alone (registry reachable) | Fix 1 alone (additive propagator) |
| --- | --- | --- |
| Scored panel rows moved | **6 of 42** — pea Hexanal ×4.3173, soy Hexanal ×9.5401, soy Nonanal ×1.0667, across the four internal snapshots | **23 of 42** — the sulfur channels rise, the pyrazine / disulfide / furfural rows fall by ×0.9747 to pay for them |
| Hofmann MFT / FFT | unchanged | 242.38 → **283.59** / 217.99 → **297.28** |
| External hold-out (8 points) | **bit-identical** | **bit-identical** |
| Matrix ln-sigma | bit-identical | bit-identical |

**Both fixes leave the external hold-out untouched, and the reason is worth stating.** All four
hold-out bundles execute the `matrix_only` path: it passes compound *names* to the registry
(so Fix 2 never applied to it) and it bypasses `predict_from_steps` entirely (so Fix 1 never
applied to it). **The external hold-out therefore tests the lipid-oxidation and observability
lane and says nothing whatsoever about the Maillard network propagator** — an eight-point
insensitivity that no previous wave had measured.

**What got worse, and stayed worse.** The flagship `cys_ribose_140C_Hofmann1998` contract now
fails on *both* of its criteria (1.4864 / 0.1267 dex against 1.45× / 0.09 dex) where Wave P had
it failing on one. FFT was already over-predicting and rose 1.36×. Nothing was refitted to
recover it: the two lanes share their upstream trunk, so a barrier that pushed FFT down would
take MFT with it. The Monte-Carlo external-literature interval **widened** 0.8495 → 0.9463 dex
while its coverage stayed at 1/3 — extra width that bought nothing.

**And one honest loss of signal.** Making the registry reachable removed the
`mechanistic_blocker` flag from the internal snapshots' Hexanal and Nonanal rows, emptying the
mechanistic-refinement watchlist (2 → 0 benchmarks). The rows now resolve a real record whose
evidence strength is `process_state_mismatch` — the registry's honest label for a factor
reached only through the `intermediate_matrix → ambient_slurry` fallback — but
`_matrix_closure_action` has **no branch for that label**, so it is scored exactly like a
genuine class-anchored transfer. The governing decision did not change (still
`hold_observable_first`, still zero approved DFT jobs), so no compute was unlocked; what was
lost is a warning. Pinned at 0 with the cause, and carried as an open item rather than repaired
in the same pass as the fix that exposed it.

### What Round 3 means for the rest of the repository

The 6 findings above came from ~24 checkable claims in the *live panel* — the most scrutinised
data in the repository, already through two audit rounds. There are roughly **200
content-unverified anchors** in `data/lit/`. On this base rate they must be presumed similarly
contaminated until someone reads their sources. That is the honest generalisation, and it is
not comfortable: **the panel is not the worst-audited part of this repository, it is the
best-audited part.**

Five benchmarks remain unverifiable without library access — Trikusuma 2020 (the last
unchecked pillar of the matrix lane, and the only fit-recovery benchmark still scoring
`pass`), Bi 2020's tables, Liu 2023, Bolton 1994 and Ramirez-Jimenez 2000.

### Wave S1b — the pH and water-activity physics that was written but never connected (2026-08-27)

Wave S2 built the repository's first **ordinal** accuracy measurement — a 52-claim
directional panel — and found the model at **16/19 on sugar identity, temperature, time and
precursor loading, and 2/10 on pH and moisture**. Worse than a coin, and systematically so.
Wave S2 also diagnosed the cause by inspection, and it was not a modelling failure: it was
three **routing** defects. Wave S1b confirmed each independently and fixed them. **No
correction curve was reshaped, no constant refitted, no barrier moved, and the panel was not
iterated against** — the wiring was chosen from the chemistry, measured once, and reported
whichever way it came out.

1. **`ReactionConditions.get_ph_multiplier` had never been executed on a prediction.** It
   encodes the 1,2- vs 2,3-enolisation route selection, it is unit-tested for the correct
   signs, and `grep` finds its only callers in `kinetics.py`, `pathway_ranker.py` and
   `cantera_export.py` — none of which is reachable from `evaluate_benchmark_payload`. It is
   now called once, from `get_rate_constant`, gated to the three enolisation branch-point
   families so it cannot double-count against the amine-ionisation term.
2. **The pyrazine branch of `_ionization_correction` keyed on the substring `"pyrazine"`.**
   Measured: of the **29 distinct reaction families this engine emits** across
   `data/benchmarks/`, not one contains that substring. The branch returned 1.0 at every pH,
   and the model moved 2,5-dimethylpyrazine the wrong way with pH against two independent
   direct measurements. The pyrazine step here is `Aminoketone_Condensation`; it is now named
   explicitly, at the branch's own unchanged pKa 6.5.
3. **`_water_activity_correction` reached 3 of those 29 families**, none on the furan/HMF
   track, and its dehydration branch keyed on `"furfural"` — equally dead. Membership is now
   decided by **measured net water stoichiometry**: water-releasing condensations and
   dehydrations get the unchanged `1.3 − aw` curve, water-consuming hydrolyses get mass
   action in `aw`, net-zero families get nothing.

**Seven dead substring keys were found in total** and are itemised in
`docs/validation/directional_accuracy_report.md` §A3, along with one that is live and was
worse than dead: `"condensation"` in `get_ph_multiplier`'s Schiff branch matches
`Aminoketone_Condensation`, so an ungated call would have handed the *pyrazine* step an
*acid*-peaked Gaussian. This is the same defect class as `Furanone_Strecker_Reduction`
(Wave P) and the Wave I offset keys; the fix is explicit family-name sets throughout.

**The result, reported both ways because only both ways are honest.** The directional
headline rose 19/29 → **20/29**; pH rose **2/7 → 4/7**; the non-pH/aw bucket *fell* 17/19 →
16/19. The two gained rows are exactly the two dimethylpyrazine-vs-pH claims defect 2 was
diagnosed against. **pH and moisture combined are still 4/10 — at or below chance — and
moisture is still 0/3.** Wave S2's recommendation to guard pH and moisture recommendations
at runtime stands unchanged.

**And it cost four benchmark rows**, each degraded and none clawed back: Hofmann
1.486× → **2.209×**, Bolton 748× → **6731×**, Cerny 2.787× → **23.406×**, Resconi 4.402× →
**5.457×**. Attribution was measured one fix at a time (by emptying the other two family sets
at runtime, no source edit): Bolton and Cerny are carried almost entirely by the water-activity
term, Hofmann by the route-selection term. **The Hofmann case must be stated precisely, and
Wave S1b's first attempt at it was wrong.** That draft called it "a genuine conflict rather
than a bug" — the model's FFT/furfural arm winning at pH 5.0 against a benchmark measuring
MFT 342 > FFT 200. **Wave S2b, the same day, showed there is no measurement on the other
side.** The 342 / 200 ppb targets were derived inside this repository from
`data/benchmarks/maillard_validation_benchmarks.md` §1.3 — an abstract-reconstructed range
table (`~0.02–0.05` and `~0.01–0.03` mol %) committed in the *same commit* that created the
benchmark file — and both values are interior points of two **overlapping** invented bands
(MFT 228–571 ppb, FFT 114–342 ppb). The ordering is midpoint selection, and the
1.45× / 0.09 dex contract is ~1.7× tighter than its own source band. The pH mechanism and the
degradation are both real; the yardstick is not. Wave S2b's recommendation — retire the
contract, demote the benchmark out of PRIMARY, mark the values `no_verifiable_source` — is an
owner decision and **no wave has executed it**; Wave S1b did not relax the contract either.

**`status_counts` are unchanged** (scale-gap 8 / pass-no-ranking 2 / pass 4) and strict-ready
is unchanged at 0/14.

**One self-inflicted bug, found and fixed inside the wave, reported because it changes how the
numbers read.** The first wiring left Wave H's `"enolisation"` substring guard in front of the
new water-activity sets, which silently excluded `Enolisation_1_2` — the furan track, the one
family the fix most needed to reach. It was caught by inspecting the factor table, not by the
score. The guard was removed from that function; the explicit sets subsume it strictly better,
and Wave H's own pin on `Enolisation_2_3_Amadori` still passes.

### Wave S2c — the tightest contract was anchored to the repo's own guess (2026-08-27, owner-approved)

**Round 3 addendum. One sentence carries this wave: THE SULFUR BRANCH NOW HAS ZERO ABSOLUTE
LITERATURE ANCHORS.**

Round 2 quarantined three fabricated sulfur benchmarks and the repo then said, in the README,
in `AUDIT.md`, in `src/barrier_constants.py`, in two refit records and in four test docstrings,
that **one** literature anchor survived: `cys_ribose_140C_Hofmann1998`, MFT 342 ppb / FFT
200 ppb, DOI `10.1021/jf9705983`. Wave S2b (a read-only literature wave) settled where those
two numbers came from. Neither Hofmann & Schieberle 1998 nor Mottram & Nobrega 2002 is the
source. They were derived **inside this repository**, from
`data/benchmarks/maillard_validation_benchmarks.md` §1.3 — an abstract-reconstructed range
table committed in `c7efbbc`, the *same commit* that created the benchmark JSON. That table's
row reads `| Ribose + Cys, pH 5 aqueous | 140 | 30 min | ~0.02–0.05 | ~0.01–0.03 |`. On the
benchmark's declared — and itself unattested — 10 mM basis with MW 114.17:

* MFT `0.0300 mol % × 0.010 M × 114.17 = 3.4251e-4 g/L = 342.5 → **342 ppb**`
* FFT `0.017321 mol %` (the exact **geometric mean** of the `~0.01–0.03` band)
  `→ 197.8 → **200 ppb**`

Confidence ~90%; the arithmetic is exact to the rounding, both values sit inside their bands,
and no external source containing 342 or 200 was found after exhausting every retrieval route.

**Why §1.3 is guesswork and not a transcription — four tells, none needing the paywalled body.**
(1) Its only bold, unhedged row (`1.4` / `0.05` mol %) is *verbatim from the abstract*, and
belongs to a **dry 180 °C / 6 min** intermediate system, not to an aqueous one; every other row
is tilde-hedged. (2) Its `Furfural + H₂S … ~0.5` cell is arithmetic on the abstract's "a 10
times higher efficiency" (10 × 0.05), yet is assigned 140 °C / 30 min — conditions the abstract
nowhere supports. (3) Its `Glucose + Cys` cell reads `~10× lower than ribose`: prose typed into
a numeric column. (4) The 140 °C / 30 min conditions are transplanted from §1.1, which is
**Mottram & Nobrega 2002** — a *headspace peak-area* study, as §1.1's own note says — licensed
by §1.1's sentence *"Fully quantitative SIDA data for MFT/FFT **from the same system**
available in Hofmann & Schieberle (1998) — see §1.4 below"*, whose cross-reference points at
Brands & van Boekel. That sentence is now deleted in place, quoted, and explained.

**What was executed.** Values marked `no_verifiable_source` and **left unedited** — the shipped
number is the evidence. The **1.45× / 0.09 dex contract retired**: it was ~1.7× tighter than
the 2.5× spread of the band its own target was interpolated from, which is the same
tolerance-too-narrow tell recorded for the Round 2 quarantines. Nothing looser was invented to
replace it; the global free-precursor default (1.5× / 0.10 dex) is inherited, is marginally
looser, is failed by more, and that inheritance is written into the file rather than left to be
discovered. `metadata.tier` demoted **`PRIMARY → REFERENCE`**, which removes strict-gate
eligibility. `thiol_addition_pentodiulose` **reverted 26.35 → 28.60** — the un-fitted Wave N
class value — and `results/validation/sulfur_barrier_refit_pentodiulose.{json,md}` **retracted**,
the same treatment Wave I gave `hydrolysate_observability_rederivation` when its two fit targets
turned out to be fabricated. §1.3 kept under a dated **ABSTRACT-RECONSTRUCTED GUESSWORK**
warning: it is the provenance record of a fabrication that reached shipped code, and deleting it
would destroy the evidence.

**What it cost, published rather than absorbed.** Removing a fit to a non-measurement makes the
model look worse against that non-measurement, and it should: Hofmann MFT **154.85 → 78.09 ppb**
(2.2086× → **4.3797×** under), FFT **267.50 → 293.67** (1.3375× → **1.4684×** over), MALE
0.2352 → **0.4041 dex**; Cerny 2008 MFT 23.406× → **25.741×**; and the pentose ≫ hexose ordering
headline **18.27× → 8.26×**. Every aggregate held: strict-ready **0/14** (a *failing* contract
was retired), honest literature coverage **0/3**, MC coverage 28/35, and all eight external
hold-out points **bit-identical for a third consecutive wave**.

**One thing got better, and it is the thing that matters.** Only **4.27×** of the pentose ≫
hexose ordering was ever structural, and that number did not move — so the share of the claim
carried by mechanism rather than by a barrier gap went from **23% to 52%**, and what is left
riding on a gap is 1.05 kcal/mol between an *estimated class value* and an unconstrained legacy
fit, where it used to be 3.30 kcal/mol between a **fitted** barrier and that legacy fit. No part
of the surviving ordering claim now traces to a number this repository invented. The one Hofmann
claim that survives is *pentose ≫ hexose* itself, because that is stated in the paper's
**abstract**, which is retrievable — unlike the yields table, which is not.

**What is still open, stated plainly.** The benchmark still declares `measured_volatiles`, so
its rows are still enumerated and still scored; the fold errors above are errors against a
fabricated yardstick and must never be quoted as accuracy against literature. Taking it out of
the scored population entirely — via `reference_volatiles` (the treatment
`thiamine_cys_xylose_145C_Cerny2008` already carries) or via quarantine — and rebuilding it from
an ILL copy of the paper in native **mol %** are owner decisions this wave did not take, because
both change panel membership. The ILL request pack is in `tasks/audit_remediation.md`
"## Wave S2b" §(f). If the paper turns out to have no aqueous ribose/cysteine row, the honest
outcome is to retire the absolute sulfur anchor entirely and say so.

## Wave U — the network's first out-of-sample test (2026-08-27)

Every accuracy claim in the sections above is about the *matrix* lane or about systems the
model was calibrated on. Wave S1 found the structural reason: **all four bundles in the
external hold-out run the `matrix_only` execution path.** They read a lipid-oxidation load
off a matrix profile and return before `Recommender.predict_from_steps` is ever called.
That is why three consecutive waves of reaction-network work — an additive flux propagator,
a pH/water-activity rewiring, and a shipped barrier revert — moved dozens of in-panel rows
and left **all eight hold-out points bit-identical**. The invariance was evidence about the
hold-out's coverage, not about the model. **The chemistry this repository is actually about
had never been scored on a system it had not already seen.**

Wave U built the missing thing: twelve content-verified free-precursor literature points,
frozen under `data/benchmarks/external_validation/maillard_path/`, every one of them
executing the reaction network. **Every value was read out of the actual paper or thesis and
is quoted verbatim in its bundle**, with the access route and retrieval date; nine candidate
sources were declined and the reason recorded for each (the full table is in
`tasks/audit_remediation.md` "## Wave U"). No conversion uses an assumed molar basis — the
one source reporting a basis-free ratio keeps it as the native unit, which is the direct
lesson of the retired 342/200 anchor.

**The predictions were frozen BEFORE any calibration wave saw these points.**
`results/validation/maillard_path_holdout_frozen_predictions.{json,md}` records the git HEAD
it was generated from, and it is un-gitignored specifically so that a later wave cannot
regenerate it after calibrating and compare it with itself. That file, not a re-run, is what
the pending rate-calibration work must be scored against.

### The baseline, stated plainly

| Measure | Result |
| --- | --- |
| Points / scored targets | 12 bundles · 22 targets · 21 quantitatively scored |
| **Median fold error** | **6.04×** |
| Worst / best | 52.59× (FFT at 130 °C) / 1.52× (MFT at 100 °C) |
| Within 10× | **12 / 21** |
| Within-point orderings | **8 / 12** pairs |
| Cross-bundle response directions | **3 / 6** |
| Structural zeroes | **1** — the model returns *nothing* for glucose with no amino-acid partner |

**6.04× is better than this repository's own in-panel numbers would have predicted, and the
encouraging part is real: the sulfur branch, which Wave S2c established has ZERO absolute
literature anchors, predicts MFT at 100 °C to 1.52×.** That is an uncalibrated branch
landing inside a factor of two on an isotope-dilution measurement it had never seen.

**And the median is the least interesting number here, because the shape is wrong in three
specific, nameable ways that a fold-error median cannot show:**

1. **The sulfur branch has the temperature dependence backwards.** Measured MFT *falls*
   4.0× from 100 °C/4 h to 130 °C/0.5 h; the model has it *rising* 4.55×. Response ratio
   18.3, direction wrong. The good 1.52× at 100 °C and the 12.1× miss at 130 °C are the same
   error seen from two ends.
2. **Acrylamide barely responds to time.** Measured 28 → 1459 ppb over 10 → 30 min at
   180 °C, a 52× rise; the model moves 61 → 76 ppb, a 1.24× rise. Response ratio **0.024** —
   the model is roughly forty times under-responsive in time on this route.
3. **Two of three pH responses point the wrong way.** From pH 5 to pH 8 the source measures
   furfural ×1.47 and HMF ×1.76; the model gives ×0.774 for both. Only furaneol (×5.11
   measured, ×2.48 predicted) moves correctly.

Also frozen: the two acrylamide lanes in this repository **disagree with each other by
roughly 480×** on the same 0.2 M glucose/asparagine system (network 1143 ppb, lumped
`predict_acrylamide` 544 870 ppb). Nothing in the panel had ever put them side by side.

**What this does not license.** Twelve points is not a validation. Five of the twelve are
acrylamide or HMF systems whose route carries a declared partial contamination — the
`safety_risk_acrylamide` barrier is a transfer of the Knol group's lumped Ea — so their
*temperature* dependence is not out of sample even though their yields are; each bundle says
so in a `contamination_disclosure` block. And the median excludes the structural zero, which
is a total miss with no finite fold error.

## Wave S3 — the first rate-level calibration, and what it proved about the model (2026-08-27)

Every barrier in this repository had, until this wave, been an endpoint fit, a literature
midpoint, or an estimate by analogy. **Nothing had ever been fitted to a measured
trajectory.** Wave S3 fitted the sugar trunk to 176 point-verified concentration-time values
from Martins' glucose/glycine multiresponse experiments at 80/100/120 °C, plus a second
experiment on *isolated* Amadori compound that separates its decomposition from its own
formation. The fit corpus is `data/lit/timeseries/` and nothing else; it is machine-declared
in `results/validation/trunk_rate_calibration_refit.json` and read by the fit-target gate.

### The credibility test passes

The fitted constants were checked against **Brands (2002)**, a genuinely independent fitted
kinetic model — different amine (protein-bound lysine, not free glycine), different
sugar:amine ratio, different author, different data — none of which entered the objective.

| step | fitted @120 °C | Brands @120 °C | agreement |
| --- | ---: | ---: | ---: |
| condensation (sugar + amine → Amadori) | 8.04e-05 L/(mmol·min) | 2.4e-04 | **3.0×** |
| total Amadori degradation | 0.189 /min | 0.2805 | **1.5×** |

Two laboratories, two amines, two independent fits, agreeing to a factor of 1.5 on the rate
the whole cascade is named for.

### A standing contradiction, resolved — and neither file was right

Since Wave I this repository has carried two barrier tables that disagreed by ~6.6e8 about
which of the first two Maillard steps is rate-determining. The Martins data measure the
Amadori intermediate directly, so the data could decide. On a pseudo-first-order basis at the
experiment's own 200 mmol/L glycine and 100 °C, the fitted ratio is **44.9 (95% interval
40–45), Amadori faster** — the condensation is rate-determining.

- `src/barrier_constants.py` asserted the ratio was ~5e-05 (**wrong sign**).
- `data/lit/arrhenius_params.yml` asserted ~3.3e+04 (right sign, **wrong by ~700×**).

Both are now annotated with the resolution and the arithmetic. **Neither file's values were
changed**, and the reason is the next section.

### The finding that matters most: the model's accuracy problem is not in its barriers

Tracing the screening lane end to end established that **it consumes barriers only as
branching ratios**. The predicted magnitude is `budget × mole fraction × MW`, and the budget
(`src/projection.py:170-223`) never sees a barrier — it is built from precursor molarity,
duration and a separate apparent activation energy. A uniform shift of *every* barrier by
+2, +5 or +10 kcal/mol changes the total predicted ppb by under 0.4%. Two consequences fell
out of the same trace: the per-family Arrhenius prefactors in `arrhenius_params.yml` **cancel
exactly** and are dead in this lane, and the one place an absolute rate survives
(`recommend.py:1216`) uses a family-agnostic default.

**So this wave predicted, in writing and before scoring, that a rate calibration could not
move the hold-out** — and then measured it. Scored against the pre-registered baseline frozen
at commit `12f43dd`:

| | as shipped | if the derived barriers were applied |
| --- | --- | --- |
| Hold-out targets moved | **0 / 22** | 21 / 22 (11 better, **10 worse**) |
| Hold-out median fold error | **6.0388× — unchanged, exactly** | 7.6110× (26% worse) |
| Directional panel | **21/29 — byte-identical** | 18/29 |

All three of Wave U's named structural errors are unmoved, as predicted: the sulfur
temperature dependence is still backwards, acrylamide is still ~40× under-responsive in time,
and furfural and HMF still move by an identical pH factor. That last one gained *positive*
evidence: under the counterfactual both shift to exactly 0.3008, still identical to each
other after a large and differential barrier change — confirming Wave U's inference that they
share a single pH multiplier.

**The honest headline is therefore a negative result with a positive cause.** The trunk's
rates are now measured, they agree with an independent laboratory to 1.5×, and applying them
to the screening lane would make the model *worse*. The absolute-accuracy deficit lives in
the projection budget and in missing chemistry — the network has no glucose→fructose
isomerisation, no formic or acetic acid, no melanoidins, and no caramelization lane — not in
the barrier table. `FAST_BARRIERS` is unchanged by this wave; the derived values are
published with their arithmetic as an owner decision.

One value did gain independent support: `amadori_rearrangement` = 23.0 kcal/mol derives to
**23.20** from the fitted rate. It is the first and only constant in that table with evidence
from a measured trajectory behind it.

## Wave S4 — binding physics as an observability mode: measured constants, zero fitted parameters (2026-08-27)

Wave S3 showed the absolute-accuracy deficit does not live in the barriers. This wave tested
the next layer down the matrix lane: the **observability factors** that convert a predicted
total into the number a paper reports. Every one of them was fitted — several back-solved
from the benchmark they are then scored against, and the 1-hexanol pair (`0.143 / 0.063`)
back-solved from two values Wave T3 proved appear in no publication.

`src/protein_binding.py` computes the same factor from measured protein-binding data
instead: `f_free = 1/(1 + a_p·Pow·c_p)`, a single-site Langmuir in the dilute-ligand limit,
with `a_p` in litres per gram of protein. **Zero fitted parameters.** Every constant is
transcribed with a verbatim quote into `data/lit/binding_constants.yml` from a full text that
was retrieved and read. It ships as a *mode*, selectable per run; `fitted_factors` remains
the default and no shipped number moved.

**The model reproduces percent-bound measurements it was never built from.** `a_p` was fitted
by one laboratory on 2-alkanones by APCI-MS; the check below is against other laboratories,
other methods and other chemical classes, with nothing adjusted:

| measurement | measured % bound | predicted from `a_p` | residual |
| --- | ---: | ---: | ---: |
| Wang 2015 thesis, pea isolate, 2-octanone, pH 7 | 31.90 | 32.64 | +0.74 pts |
| Wang 2015 thesis, pea isolate, 2-heptanone | 13.90 | 16.49 | +2.59 pts |
| Heng 2005 thesis, pea vicilin, 2-heptanone | 19.00 | 19.27 | +0.27 pts |
| Heng 2005 thesis, pea vicilin, 2-octanone | 33.00 | 36.95 | +3.95 pts |
| Barallat-Pérez 2023, plant isolates, octanal | 52.76 | 55.61 | +2.85 pts |
| Wang 2015 thesis, pea isolate, octanal, pH 8 | 61.87 | 50.88 | −11.0 pts |

**And then it does not rescue the matrix lane, for a reason worth stating.** Scored on the
eight-point external hold-out — never fitted, by anyone — against the incumbent and against a
*null model* that applies no observability factor at all:

| mode | hold-out median fold | CI coverage | in-panel median fold |
| --- | ---: | ---: | ---: |
| `fitted_factors` (shipped) | 93.68× | 3/8 | **1.0004×** |
| `unit_observability` (null) | **67.42×** | 4/8 | 5.92× |
| `binding_physics` | 68.18× | **5/8** | 5.92× |

Read the two comparisons separately, because they say different things. **The fitted factors
beat everything in-panel and lose to doing nothing out of sample** — the textbook signature of
constants fitted to the rows they are scored on. **The binding physics is a wash**: it moves
the median by less than 1% against the null model, because on 12 of 16 hold-out rows there is
no usable binding datum and the mode reduces to the null model exactly. On the two scored rows
where it genuinely applies it is 1–1: Liu nonanal improves 94.2× → 1.85×, Liu hexanal degrades
48.2× → 121×.

**The sharpest single result is a unit argument, not a score.** An observability factor is a
fraction of a total and cannot exceed 1. The shipped constants are **4.32 (pea) and 9.54
(soy)**. Whatever those numbers are, they are not observability — they are absorbing an
absolute-scale deficit that lives in the marker yields, which were themselves built from the
retired 260 / 380 ppb values. No binding model can repair that from the observability side.

The Liu-vs-Pratap-Singh question Wave R left open is now answered quantitatively and the
answer is negative: binding predicts Liu's water-calibrated numbers under-read the true total
by **2.51×**, which widens the Liu-vs-Pratap gap from 9.95× to **24.9×** instead of closing
it. The two are a materials difference, not a method artefact.

Full row-by-row result: [matrix_binding_mode_comparison.md](results/validation/matrix_binding_mode_comparison.md).

## Wave S5 — the front door made to match the evidence, and nine documents fewer (2026-08-28)

Every wave before this one improved what the repository *knows*. This one changed what it
*says first*, because the two had drifted: the tools led with absolute ppb numbers, and every
out-of-sample measurement here says those are wrong by 6x to 94x, while the ordinal claims are
right 21 times in 29.

**One entry point, ratios first.** `python scripts/maillard.py` has three verbs — `compare`
(per-compound A/B ratios, dominant pathway, per-axis reliability), `predict` (a range, not a
point, with the caveats inline) and `rank-experiments` (the value-of-information queue). It is
orchestration over existing machinery; no science changed to build it. Absolute numbers appear
only behind `--absolute`, and their caveat prints in the same block. The per-axis reliability
tags are read from `docs/validation/directional_accuracy_report.md` **at runtime**, so a
sugar swap reports `trust (8/8)` and a pH sweep reports `do-not-use (4/7)` — and re-scoring the
panel changes what the tool tells a user with no code change.

**The model card is generated, not written.** `scripts/generators/generate_model_card.py` reads
the directional panel, both hold-out artifacts, the mode comparison, a live re-run of the
benchmark panel, a recount of the `no_verifiable_source` census, the sulfur benchmark's anchor
status and the three CI gates, and writes a claim-type x system-class validity domain into
README.md between markers. `--check` turns drift into a failing command.

**The reason it is generated is a defect this wave found.** The directional report was still
publishing **20/29** as its headline; the tree has scored **21/29** since Wave T4 moved SUG-12.
Three waves published a number one point below the truth in the file whose whole purpose is to
state that number. A second, subtler version was caught while wiring the tags: the report's
category table and its bucket table have different denominators, so summing categories to get
"excluding pH and water activity" gives 23/25 where the headline-comparable figure is
**17/19**. Both tables are now labelled and the arithmetic is written out in place.

**Nine documents were folded into their nearest living home and deleted: 30 -> 17 documents,
9157 -> 7387 lines.** Six were byte-identical duplicates of files under
`data/Gemini_Deep_Research/` (`cmp`-verified), one of them sitting under `docs/` with no warning
header at all. `docs/architecture.md` and `docs/reference/SMIRKS_SYSTEM.md` were folded into
README **with their pre-audit trust claims removed rather than relocated** — the first opened
with "High trust — use freely", the second called two templates "Validated vs. benchmarks" at
0/14 strict-ready. Eight citation strings were repointed at the surviving copies, two of which
had been dangling before this wave. `pathways.md`, which exists nowhere and is cited by three
runtime-critical data files, was **not** repointed: it is named and flagged in place, because
substituting a plausible source for an absent one is the defect this remediation exists to
remove.

All three gates pass. The citation gate's file count falls 98 -> 92 and its DOI count 1013 ->
909, entirely from the deleted duplicates; the same DOIs remain under
`data/Gemini_Deep_Research/`.

## Wave Y — the budget workstream: three findings, one fix, and two falsified predictions (2026-08-28)

Three earlier waves triangulated the same address for the model's absolute-accuracy deficit.
Wave S3 measured that the screening lane consumes barriers only as branching ratios. Wave S4
made a unit argument: an observability factor is a fraction of a total and cannot exceed 1, and
the shipped ambient hexanal factors were 4.32 and 9.54. Wave X measured 0.518 dex on a single
step against 0.952 dex end to end. This wave maps the layer all three point at, moves the one
constant the argument licenses, and reports what the move did not fix.

### The layer, measured

    observable_ppb = L(pool, lane, T, t) x Y(compound) x release(compound, lane) x cal(compound, lane, state)

`Y` and `cal` are **perfectly degenerate on any single lane** — only their product is
identified — so which one carries a scale is a *convention*, and the convention the registry
declares is that the pea ambient slurry is the reference lane whose factors are 1.0. Wave O
wrote its shared scale into `cal`, which broke that convention and produced a fraction-observed
of 4.32. There is a second degeneracy that decides the direction of the fix: `Y` also multiplies
`hydroperoxide_scale = 1.0e6`, an arbitrary constant, so **a yield above 1 is not a unit
violation and a fraction-observed above 1 is.**

### Old → new, with provenance

| constant | before | after | kind |
|---|---:|---:|---|
| `MATRIX_BENCHMARK_BASE_MARKER_YIELDS['Hexanal']` | 0.205 | **0.885036** | Wave O's shared scale, relocated. One parameter, two verified anchors, the same 1.0113× residual. |
| `pea_iso / ambient_slurry / hexanal` | 4.31725 | **1.0** | reference lane restored; independently the value Wave S4 (c) evidenced from Pratap-Singh's matrix-matched methods |
| `pea_iso / heated_matrix / hexanal` | 0.228776 | **0.0529912** | propagation of the Trikusuma back-solve |
| `soy_iso / ambient_slurry / hexanal` | 9.54007 | **0.453 / 0.205** | the soy-vs-pea ratio, restored as a self-documenting expression |
| `soy_iso / heated_matrix / hexanal` | 2.80478 | **(0.453/0.205)·(1−0.7060)** | composition rule, unchanged, following its baseline |
| `2-Pentylfuran` yield 0.502 | — | **not moved** | re-derives to 0.5017897, 0.000182 dex away — below the 0.01 dex materiality floor |
| `1-Hexanol` / `Nonanal` yields | — | **not moved** | no anchor of any kind; not fitted |
| `reference_conversion_time_min` | 12589 | **not moved** | see below |

**Nothing predicted moved on the calibrated tier.** The eight hold-out points are unchanged to
six significant figures (four hexanal rows drift 3e-7 to 8e-7 from rounding, four non-hexanal
rows are bit-identical), and the pea lane still reads 1125.278 ppb.

### What DID move, and it is the point

`_uncalibrated_prediction_ppb` reads the yield and never reads an observability factor. Wave O
recorded that no refit of *those* constants could ever move the matrix sigma derivation; a
**yield** refit is the case that sentence does not cover, and this is the first one. Soy ambient
hexanal goes from **9.4334× under to 2.1851× under** in that tier, and the scatter narrows
(centred sd 3.0951 → 2.8452 in ln space) — but the **bias grows** (3.31× → 5.95× over) because
the same relocation makes the heated lanes over-predict by 4.32× in a tier that has no
observability to correct them with. Net RMS ln-sigma 3.0166 → 3.1075; shipped 2.86 still inside
the re-derived 90% interval [2.088, 6.494]. **The prior was not re-fitted** — doing so in the
wave that moved the residuals would make the two unattributable.

### The Wave S4 claim: PARTIALLY CONFIRMED

Two factors came back under 1 (pea ambient → 1.0, soy heated → 0.6497). **Six did not, and
every one of them is a soy factor.** The reason is structural: a marker yield is shared across
matrices, so it can absorb a *global* scale error and never a *lane* one. Pinning observability
to the 1.0 that Wave S4's own evidence requires and fitting two yields against four verified
rows leaves an RMS residual of **0.3103 dex**, entirely soy-vs-pea — required-yield ratio
**2.1606× on hexanal, 5.9221× on 2-pentylfuran**. Those two differ from each other, so the
deficit is *compound*-specific as well as *lane*-specific and **cannot be repaired by the soy
lipid profile either** (that would move both linoleate markers by one factor). Wave S4's
diagnosis was correct about the global scale and incomplete about the rest.

### The projection budget: a hard new constraint, and a decision not to apply it

Wave X's single-step systems print molar yields on the fed precursor, which is exactly what the
budget's `conversion_extent` predicts — with no allocation term in between. That gives the
repository its first constraint on the budget from *below*:

> `hofmann1998_furan2aldehyde_h2s` measures **11016 ppb** of FFT. The model's entire volatile
> budget at those conditions, expressed in that analyte's own molar mass, is **2413.65 ppb**,
> and the allocation already gives that analyte **100 %** of it. The row is unreachable at any
> allocation, any barrier and any observability factor. Required budget scale **4.564×**.

Three rows are unreachable in this sense, all of them Wave X step-level rows. And yet the same
table demands the opposite elsewhere: required scale spans **0.0022× to 4.564×, a factor of
2075**. Grouped by how many steps the route actually took, short paths (≤2 steps, 5 rows) want
the budget **1.07×** and long ones (≥4 steps, 11 rows) want it **0.0185×** — a **57.9×** split.
Applying the binding constraint costs the panel (mean |log10| 0.9283 → 1.0119 dex, 8 rows
better and 8 worse), so `tau_ref` is **not** moved. **The deficit is not a scale. It is that the
budget is allocated from the limiting precursor with no regard for how many steps separate that
precursor from the product** — which is Wave X's step-vs-cascade result with a mechanism under
it. Record:
[projection_budget_step_yield_constraint.md](results/validation/projection_budget_step_yield_constraint.md).

### Two pre-registered predictions were falsified, and both falsifications are findings

Predictions were written to
[wave_y_prereg_predictions.json](results/validation/wave_y_prereg_predictions.json) before any
constant was edited.

* **Y-P1 — "0 of 22 maillard_path hold-out targets move" — FALSIFIED: 8 moved.** Not by this
  wave. Ablation: with Wave Y's five constants reverted in memory, the same 8 still differ;
  with Wave X's norfuraneol channel replaced by a no-op, **all 22 are bit-identical to the
  pre-registration**. Wave X priced that channel against the Wave W panel (0.9241 → 0.9518 dex)
  and never re-scored the frozen hold-out; the unmeasured price is **8/22 targets and a median
  fold error of 6.0388× → 10.8638×** on the repository's only free-precursor out-of-sample
  surface. Wave Y's own contribution is **0 targets**.
* **Y-P8 — "the four synthetic snapshots are bit-identical" — FALSIFIED: every network-derived
  volatile on them moved by exactly ×1.27604.** The mechanism is a genuine cross-layer leak: on
  the `matrix_precursor_augmented` lane the marker yields set the *injected* marker molarity,
  which enters `_relative_precursor_load_factor` — a geometric mean over all positive initial
  concentrations — and therefore the volatile budget. 4.317249^(1/6) = 1.27612 with six
  positive species. **A lipid-oxidation observability constant is setting part of the Maillard
  volatile budget.** Carried as an open item; it is a structural change, not a constant.

Everything else scored as predicted: the Wave X step-level rows are unchanged at 0.5177 dex, the
hold-out is unchanged, the honest external-literature coverage is unchanged at 4 hits, and the
matrix sigma moved exactly as predicted in mechanism (both uncalibrated hexanal predictions
×4.317249) though not in aggregate direction.

### What Wave Y cost

| aggregate | before | after | reading |
|---|---|---|---|
| Panel passes, strict | 4/14 | **0/14** | the last four were synthetic snapshots; ×1.936 inherited from Waves W/X, ×1.276 from this wave, against a ×2.00 threshold |
| Panel passes, lenient (`benchmark_summary.md`) | 6/14 | **2/14** | the two Pratap-Singh `pass-no-ranking` fit-recovery rows survive |
| MC panel 90 % CI coverage | 32/47 | **28/47** | all four lost hits are synthetic rows |
| Honest external-literature coverage | 4 hits | **4 hits** | unchanged — the headline that carries evidence did not move |
| Matrix hold-out median fold | 93.68× | **93.68×** | unchanged |
| Wave X step-level mean \|log10\| | 0.518 | **0.5177** | unchanged |

The snapshots were **not** regenerated. Refreshing them would restore 0/14 → 4/14 and in the
same motion absorb the Waves W/X drift that nobody has reported.

### The barrier sensitivity deficit, decomposed (on a Wave Z1 lead)

Wave Z1 measured a ~13× shortfall between the Arrhenius expectation of a 2.30 kcal/mol barrier
move (15.9×) and its effect on the Hofmann Table 4 objective (1.22×). Reproduced here to
1.2179×, and split into **two** stages rather than one:

| probe | analyte share of budget | Arrhenius | flux ratio | ppb ratio | stage 1 (span) | stage 2 (allocation) |
|---|---:|---:|---:|---:|---:|---:|
| `furan2aldehyde_h2s` / FFT | 0.992 | 15.93× | 1.4121× | **1.0000×** | 11.28× | 1.41× |
| `norfuraneol_h2s` / MFT | 0.465 | 15.93× | 6.5817× | **1.2179×** | 2.42× | 5.40× |
| `cys_ribose_140C` / MFT | 0.109 | 16.47× | 1.7502× | **1.5982×** | 9.41× | 1.10× |

Stage 1 is the log-sum-exp **series span**: the largest barrier on a route dominates it, so
moving a smaller one barely moves the flux. Stage 2 is the **allocation normalisation**: a
channel that already holds most of a fixed budget cannot grow. On the first probe a 15.9× rate
change produces *exactly zero* change in predicted ppb, and the sum of predicted ppb is
invariant to six figures. **Barriers set the split; the budget sets the scale.** And the budget
work cannot fix it — under the counterfactual `tau_ref` every ratio in that table is identical
to 3e-16, because a common factor cancels out of a ratio. Said plainly: *the normalisation is
deeper than the budget constants.* Record:
[barrier_sensitivity_deficit.md](results/validation/barrier_sensitivity_deficit.md).

### One decision deliberately NOT taken

The Trikusuma heated-pea nonanal factor is stale by exactly `linoleic/oleic = 50/22 = 2.2727`,
because Wave P moved nonanal to the oleate pool without propagating the constant. Wave Y
re-derived that identity independently from the layer map and then found Wave P's own comment
block recording the same diagnosis and an explicit refusal to refit — on the grounds that
re-solving it would re-absorb the correction that exposed it, and that Trikusuma is still
content-unverified. Both grounds hold. **The refusal stands; Wave Y contributes only the
independent confirmation.**

## Wave R1 — the barrier table the model was not using (2026-08-28)

**The finding, in one sentence: the repository contained an automatic search that tuned
reaction barriers to improve the benchmark panel the model is then scored against, wrote
the result into a tracked data file that `get_barrier()` applies to every prediction, and
was invoked from the middle of the repository's own documented regeneration command.**

This is the same defect class as every other finding in this document — fit, then score,
without declaring the fit — but it is the first one found *in the code that produces the
barriers themselves*, and the first one whose arming mechanism was a `make`-style command
the project tells people to run.

### The mechanism

`src/barrier_constants.get_barrier()` reads `accepted_offsets` from
`data/lit/refinement_surrogate_patches.json` and **adds** the matching entry to the audited
`FAST_BARRIERS` value. `src/refinement_campaign.build_refinement_impact_artifact()` filled
that map by trying candidate offsets, keeping any that lowered the panel's total score, and
returning them as `patch`. `scripts/generators/generate_refinement_governance.py` ends with
`patch_path.write_text(json.dumps(impact_payload["patch"]))` — a wholesale overwrite of the
tracked file — and that generator is line 181 of `scripts/docker_maillard.sh`, inside
`scientific_lane()`, **ahead of** the figure generators, the validated-envelope report and
`scientific_fast_lane()`, which is the scientific-regression pytest run.

Three properties make it worse than a stale constant:

1. **It is circular.** The acceptance criterion is the score on the evaluation panel. The
   panel is then published as evidence. Nothing declared this to
   `scripts/ci/fit_target_gate.py`, so both gates passed throughout.
2. **It is saturated.** All nine accepted offsets were exactly ±3.0 kcal/mol — the search
   bound. A fit pinned to its bounds reports the bound, not an optimum. The true optimum,
   if one exists, is outside the search space and was never looked for.
3. **It is large.** At 150 °C, 3.0 kcal/mol is a ~35× rate factor.

| family | armed | audited table | error |
| --- | ---: | ---: | ---: |
| `Schiff_Base_Formation` | 18.0 | **15.0** | +3.0 |
| `Retro_Aldol_Fragmentation` | 29.0 | **32.0** | −3.0 |
| `Thiol_Addition` | 31.6 | **28.6** | +3.0 |

`thiol_addition_pentodiulose` — the constant Waves N/P/S2c argued over at length — was
**not** reached by any offset and shipped at 28.60 throughout. The `thiol_addition` *class*
value is the one that was silently 31.6.

### What was contaminated, and what was not — measured, not assumed

The tracked patch file at git `HEAD` carried an **empty** map. So a clean checkout shipped
the audited table, and every artifact in this tree was produced in that state. Re-running
the full generator sequence after the retirement reproduces `benchmark_summary`,
`external_validation_report`, `prediction_uncertainty`, `matrix_sigma_residual_derivation`,
`experiment_value_ranking`, `loo_leverage`, the family overviews and the figures
**byte-identically**; the free-precursor hold-out re-run reproduces Wave Y's shipped
numbers to the last digit. **The published accuracy numbers were not inflated.**

What was live is the *arming*: anyone who ran `scientific_lane()` — the repository's own
instruction — silently replaced three barriers and then ran the scientific regression suite
against the replaced model, with no message and no artifact recording it. Whether a given
local tree was honest depended on whether someone had checked the file out since. Running
the pre-retirement code today re-accepts the same nine offsets, so this was a live trap,
not a historical one.

### The price, on the pre-registration

The `maillard_path` hold-out (frozen at `12f43dd`, never available to any fit) scored under
both states — `results/validation/holdout_prepost_barrier_offset_retirement.md`:

| metric (32 targets) | fit armed | fit retired (ships) | direction |
| --- | ---: | ---: | --- |
| median fold error | 10.05× | **10.86×** | worse |
| median \|log₁₀\| | 1.0020 dex | **1.0360 dex** | worse |
| worst fold error | 565.2× | **506.4×** | better |
| within 10× | 15/31 | 15/31 | unchanged |
| targets differing from the 12f43dd pre-registration | **21 of 22** | **8 of 22** | the result |

**The honest model is ~8 % worse on the median of a hold-out the fit never saw.** That is
the expected sign: a panel fit that generalises slightly is still a panel fit. The last row
is what settles it — armed, the model contradicted its own pre-registration on 21 of 22
targets; retired, on 8, and those 8 are Wave X's norfuraneol channel, already attributed.

The eight-point **matrix** hold-out moved **0 of 8 points**, because it runs `matrix_only`
and never reaches the reaction network. A hold-out that cannot detect a 35× rate error in
three barriers is not evidence about barriers.

### What was done

* `accepted_offsets` permanently empty, with a retirement note recording all nine retired
  values. The note is emitted from `src.refinement_campaign.RETIREMENT_NOTE` rather than
  hand-written into the JSON, because the generator overwrites that file wholesale and a
  hand-written explanation would survive exactly until the next regeneration.
* The auto-acceptance path records candidates as `candidate_offsets_not_applied`
  diagnostics — the search is still a useful attention pointer — and applies nothing.
* `tests/unit/test_wave_r1_barrier_offset_retirement.py`: 57 test cases, including
  `get_barrier(family) == FAST_BARRIERS[family]` for **every one of the 51 families**,
  because the retired offsets reached families by three different routes (exact key,
  normalised alias, and the `schiff`/`amadori`/`enol`/`strecker`/`cys` substring map inside
  `get_barrier`) and a sampled guard would have to guess the next route.

### What this generalises to

Every previous round of this audit looked for fits in *data*. This one was in a
*generator*, and the generator was in the *build*. The lesson for the remaining `[P]` list
is that "does any script write into a tracked file that a runtime path reads?" is a
question worth asking of every generator in `scripts/generators/`, not only of the ones
with `refit` in the name.

## What this repository is now for

A research prototype of the *architecture* alternative-protein flavor science needs —
reaction enumeration with provenance tiers, matrix-aware observability, uncertainty
propagation, and a value-of-information loop that converts every model failure into a
ranked wet-lab request — together with a documented map of the kinetic measurements the
field does not have. It is not a validated quantitative predictor, and its reports say so
on every surface.

## Cleaning branch, 2026-09-01/02 — what moved and why (dated addendum)

Numbers above are as of their own dates. The cleaning branch (Phases 0–4a, plus the owner's
round-2 decisions) changed the following, each as a deletion or a relabelling, never a
verification: the QM/DFT lane and `data/qm/` (census 120 → 102); the mocked
`protein_source_registry.json` withdrawn with its code (102 → 87; no default-run output
moved); the two `*_ProtocolPilot2026` benchmark files, byte-identical twins of the
`*_Internal2026` synthetic snapshots, deleted together with the hexanal/nonanal "closure"
lane that compared the twins (panel 23 → 21; MC panel 20 → 18 benchmarks, 47 → 35 rows;
`honest_literature_coverage` 4/13 unchanged; evidence-role split 14/5/4 → 14/5/2);
`uncertainty_posture` / `validated_status` labels dropped (read by nothing); the two
literature ledgers moved to `results/literature/`; `data/Gemini_Deep_Research` renamed
`data/research_corpus`. The record of each step is `tasks/data_restructure_plan.md`.

## Open items

The ledger's `[P]` entries are the live list. Two of them are stated in full here because
each one is a *retraction of something this document previously asserted*, and a reader who
stops at the summary tables would otherwise carry away the withdrawn version.

### 1. Wave H's sulfur conclusion is withdrawn, and the constant it produced is orphaned

**What is withdrawn.** Wave H's headline diagnosis — that the sulfur residual is a
volatile-budget *allocation* deficit that **no barrier value in any defensible range** can
fix — is invalid. It was measured against a network in which the MFT route was throttled by
a defect: the norfuraneol + H₂S step drew a reducing equivalent from a pool whose only source,
in a ribose/cysteine system, was pyrazine chemistry, so MFT was structurally starved
*regardless of barrier value*. A refit run against a starved route sees a flat objective and
correctly reports that no barrier helps; it cannot see why. **The saturation was an artifact
of the defect, not a property of the chemistry.**

**What replaced it, and it is a genuine win.** Sourcing that reducing equivalent from the
thiol redox couple instead (2 cysteine → cystine + 2[H], exactly atom-balanced) moved
`cys_ribose_140C_Hofmann1998` **with no barrier constant changed at all**:

| Compound | Before | After | "Measured" |
| --- | --- | --- | --- |
| 2-methyl-3-furanthiol (MFT) | 61.25 ppb — **5.58× under** | 235.32 ppb — **1.45× under** | 342 ppb |
| 2-furfurylthiol (FFT) | 61.44 ppb — **3.26× under** | 219.96 ppb — **1.10× over** | 200 ppb |

> **Correction, 2026-08-27 (Wave S2c).** The column headed *Measured* above is **not** measured,
> which is why it now carries quotation marks. Wave S2b traced 342 and 200 ppb to
> `data/benchmarks/maillard_validation_benchmarks.md` §1.3 — an abstract-reconstructed range
> table committed in the same commit as the benchmark JSON — with both values interior points of
> two invented, overlapping mol % bands. The **structural** finding in this section stands
> untouched: sourcing the reducing equivalent from the thiol redox couple really did move the
> sulfur branch by two orders of magnitude with no constant changed, and that movement is real
> whatever it is scored against. What does **not** stand is reading "1.45× under" as agreement
> with literature. The orphaned constant this section is about, `thiol_addition_norfuraneol`
> = 26.85, was fitted through a contradicted route **and** against a non-measurement; it is
> annotated rather than reverted only because no step emits its family. Its live sibling
> `thiol_addition_pentodiulose` **was** reverted, 26.35 → 28.60. See "Wave S2c" above.
> **The sulfur branch has zero absolute literature anchors.**
>
> **EXTENDED 2026-08-28 (Wave X): the panel now also holds the repository's first STEP-LEVEL
> constraints.** Six rows from Tables 3, 4, 8 and 10 of the same paper — `precursor (1 mmol) +
> H₂S or a partner (1 mmol) in 50 mL, pH 5.0, 145 °C / 20 min` — constrain *individual reaction
> steps* rather than end-to-end lumps. Five are scored (the sixth is a fit target) at a mean
> **0.518 dex**, roughly half the end-to-end panel's 0.95 dex. Nine further verified rows are
> **not executable** and are held, with the missing step named for each, in
> `data/benchmarks/step_level_unreachable/`. The wave's own result is a **rejected fit**: the
> norfuraneol → MFT step was re-added as a slow parallel channel and its barrier fit was
> refused by a pre-registered isotope gate, because Hofmann Table 4 and Cerny & Davidek's
> labelling experiment admit **disjoint** barrier ranges (≤ 23.30 vs ≥ 27.80 kcal/mol). See
> `tasks/audit_remediation.md` "## Wave X".
>
> **CORRECTED 2026-08-28 (Wave W): THIS IS NO LONGER TRUE.** The full text of Hofmann & Schieberle 1998 (`10.1021/jf9705983`) arrived by interlibrary loan and the sulfur branch now has **three** absolute, stable-isotope-dilution literature anchors — `hofmann1998_{ribose,glucose,fructose}_cysteine_145C_20min_pH5`, the pH-5.0 aqueous rows of the paper's own Table 1 (ribose FFT 121 / MFT 198 ppb; glucose 28 / 19; fructose 32 / 25, all µg per 100 mL × 10 with the volume printed in the table footnote). The paper also confirms Wave S2b's forensic finding from the primary source: 342 and 200 ppb appear nowhere in it. **The model fails all three anchors** — 12.27×, 29.58× and 14.46× worst-ratio, mean 0.92 dex — so the branch went from *unanchored* to *anchored and measurably wrong*, which is the direction of travel this audit wanted.

That this came from a structural correction and not from a fit is the entire point: it is the
only movement in this audit that is not a knob turning.

**What happened to the orphaned constant (Round 3 update).** `thiol_addition_norfuraneol` had
been left shipping at **26.85 kcal/mol**, a value fitted by Wave H against the throttled
network. Round 3 resolved that in the strongest available way: the *step itself* was retired,
because the isotope literature contradicts the route it belonged to. The 26.85 is now a
provenance-only entry — no step emits that family. Its replacement,
`thiol_addition_pentodiulose`, ships at **28.60**, the un-fitted sulfur-addition class value,
and is labelled ESTIMATED and UNCONSTRAINED rather than inheriting a number fitted through a
route that no longer exists. The Hofmann agreement fell from 1.45× to 2.25× under as a
result, and was not clawed back.

**The refit record is now stale, and was deliberately NOT re-run.**
`results/validation/sulfur_barrier_refit_hofmann.{json,md}` profiles
`thiol_addition_norfuraneol` — a family that no longer fires — and its recorded further move
(`furanone_cyclisation` 28.00 → 26.65, `thiol_dehydration` 26.80 → 27.70, objective
**0.1018 dex → 0.0558 dex**) was measured on the retired route. **It was not applied and the
generator was not re-run**, for the same reason as before and one more: stacking a refit on a
chemistry change destroys attribution, and refitting against a benchmark whose own
mol%→ppb conversion is undocumented (Round 3, finding 6) would be fitting to an unverified
target.

**Owner decision required**, and the options are not equivalent: re-run the refit against the
corrected route (better objective, permanently ambiguous attribution, and it would put the
one remaining sulfur anchor inside the fit-target set); leave 28.60 in place (honest about
provenance, but the constant stays unjustified); or re-anchor the sulfur branch on new
measurements — the only option that actually resolves it, and the one the wet-lab campaign
exists to enable. Recovering the Hofmann full text and writing down its mol%→ppb conversion
is a precondition for the first two.

### 2. The routing root cause is still live in `src/barrier_constants.py`

Wave-scale symptom fixes were applied at each *consuming* site. The defect itself was not.
In `get_barrier()`, the calibration-offset lookup runs against
`_normalize_family_key(reaction_family)`, and only *afterwards* — at the `if fm in
FAST_BARRIERS` fallthrough — is the key re-resolved through `_canonical_fast_family()`. So a
caller whose family name needs canonicalisation gets **the right barrier with the offset
silently dropped**. Verified 2026-08-27, reproducer:

```text
BARRIER_OFFSETS={"schiff_condensation": 5.0}
  get_barrier("Schiff_Condensation")   -> 20.0   # normalizes directly; offset applied
  get_barrier("Schiff_Base_Formation") -> 15.0   # canonicalizes; offset SILENTLY LOST
```

Both names denote the same reaction and both return a barrier from the same
`FAST_BARRIERS` entry. There is no error, no warning, and no test that would catch it: the
value returned is a perfectly plausible barrier. Every currently known consumer was patched
at the producing side, so no shipped number is wrong today — which is precisely what makes
this dangerous. **The next caller to use a short or aliased key re-opens the bug in silence**,
and it will present as a calibration that mysteriously has no effect. The fix is a
one-line reordering (canonicalise first, then look up the offset); it is left to the owner
because it changes the offset resolution of any alias currently relying on the broken path,
and that blast radius needs measuring rather than assuming.

### 3. `data/qm/` was never reviewed by anyone, because git could not see it

`.gitignore` line 33 (`data/*`) excluded `data/qm/` from the repository outright. The
consequence was not merely that a clean checkout failed — it was that **the directory was
outside the scope of every audit pass in this document**, including the 2026-08-26 CrossRef
sweep of 225 anchors. It contains nine "literature" barrier windows and eighteen columns of
claimed wB97M-V / revDSD-PBEP86-D4 / xTB results, and between them they carried **no DOI, no
author-year string, no input deck, no output log and no run date**. Silence read as
provenance for four months.

All 27 values are now labelled `source_status: no_verifiable_source`, which is what moved the
repo-wide count from 84 to 102 (see above). It has since moved again, 102 → 120, for the same
kind of reason: Wave T3 (2026-08-27) labelled 18 further already-shipping records — the 14
mocked protein-source profiles plus the file record in
`data/lit/protein_source_registry.json`, the two pH-release `runtime_surrogate` blocks whose
`log_slope` is an exact back-solve from an invented band, and
`ref41_ppi_sulfur_volatile_binding_v1`, which cited a reference *number* inside an LLM
research dump. **No value was added, changed or invented.** Unlike the `data/qm/` 18, all
eighteen of these are in `data/lit/` and all carry numbers, so the runtime-consumed figure
rises with them, 62 → 80 — which is to say 80 was always the true figure and 62 was an
undercount of numbers the model was already eating. Two signatures, recorded in the files themselves,
suggest the numbers were written rather than computed: every xTB and wB97M-V value falls
*inside* its own literature window, and the revDSD-minus-wB97M-V difference is non-negative
for all nine families, a clean multiple of 0.1, and confined to a 0.6 kcal/mol band. Real
functionals do not agree that tidily or that one-sidedly across nine reaction classes.

**Nothing here reaches the model** — verified 2026-08-27: the only importer of
`src/authority_benchmark_data.py` is a loader test, and `tests/benchmarks/_lane_policy.py`
merely gates the Phase 3 authority lane on the fixtures' existence. That bounds the damage;
it does not make the numbers real. Re-running or withdrawing them is an owner item.

A second-order defect worth naming: `_lane_policy.py` **skips** when the fixtures are absent.
So for as long as the data was invisible, the authority lane reported skips rather than
failures — absence of the dataset was indistinguishable from success. `scripts/ci/citation_gate.py`
now scans `data/qm/` and fails on any benchmark entry that names no source and does not admit
it has none.

### The rest

The rest of the list: the volatile-budget *allocation* layer (its
legacy heuristics are deliberately left dead rather than refitted),
re-anchoring the sulfur branch on new measurements, the **120 records marked
`no_verifiable_source`** whose numbers are now honestly unanchored — measured directly
across every tracked JSON and YAML file including nested records, of which **98 carry
numeric payloads** and **80 of those are consumed at runtime** — and the wet-lab campaign
([PPI/SPI protocol](docs/protocols/PPI_SPI_PRIMARY_BENCHMARK_PROTOCOL.md)) that would give
this model its first real calibration surface.
