# The Audit

*2026-08-26/27. This document is the narrative record of a two-day adversarial audit and
remediation of this repository. The machine-readable ledger is
[tasks/audit_remediation.md](tasks/audit_remediation.md); the evidence artifacts it cites
are tracked under [results/validation/](results/validation/).*

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
| Panel | **14 benchmarks** (2 more quarantined in Round 2); the MC uncertainty panel covers **11** of them, **35 matched rows** |
| Literature rows inside 90% CI | **1/3 evaluable** (median CI width 0.86 dex), with fitted rows removed from numerator *and* denominator |
| Fitted rows (constants back-solved from the benchmark) | 2/2 — **not evidence**; both would previously have been counted as literature hits |
| Internal-synthetic rows (reproducibility only) | 18/18 — carries zero validation weight |
| Benchmarks without blocking gaps | **0/6 predictive** (+ 1/4 fit-recovery, 4/4 synthetic; the 5/14 aggregate is all non-evidence). Fit-recovery fell 3/4 → 1/4 in Round 3 when the two Pratap-Singh benchmarks were corrected against the paper and stopped recovering. Wave O then refitted their constants onto the verified values and they recover again — at `pass-no-ranking`, which is not `pass`, so **none of these counts moved**. A refit changes the size of a recovery, never its evidential status. |
| External hold-out | **1/5 genuine extrapolations** at the pre-widening prior — and 1/5 under the wider one too since Wave O; median fold error **42.62×**, worst **2474×**. Only 4 of the 8 points are measurements at all. The median was 32.79× until Round 3 corrected two wrong table rows in the Li 2026 bundle (→ 15.31×, predictions unmoved), then rose to 42.62× when Wave O refitted the ambient hexanal observability onto the paper's verified anchor. **A refit to verified values made the out-of-sample number worse**; see Round 3, Wave O. |
| Strict-ready benchmarks | **0/14** |
| Reaction-chemistry lanes with generative templates | 5/16 (derived by enumeration, test-pinned) |
| Test suite | **1242 passed, 1 skipped, 2 xfailed, 0 failed** (2026-08-27, after Round 3). The dvipng failure carried as "environmental" for the whole audit is **gone**: `--report` now degrades to mathtext instead of raising |
| Citation gate | **0 dead DOIs as of the 2026-08-26 and 2026-08-27 sweeps; 0 waivers.** The blocking gate is **structural and offline** — DOI grammar, confabulation signatures, status coherence, repair-record completeness. It cannot detect a live DOI that resolves to the wrong paper, which is exactly how the two PMC9905368 benchmarks survived every previous check. Liveness is a **scheduled weekly sample**, advisory, not part of the required set. |

The coverage numbers *fell* during remediation — twice — and each fall is documented at
the point it happened (see the README's calibration section and the ledger). That is the
audit working as intended: each drop is the measured size of a fabricated support that was
removed. What survives at quantitative fidelity is narrow; what survives robustly is
**ordering** — pentose ≫ hexose sulfur chemistry (**3.39×** against a 3× floor; published as
15.8×, then corrected to **8.98**× on 2026-08-27, then to 3.39× the same day when the MFT
route was corrected on isotope evidence; and of that, only ~1.13× is structural — zeroing the
1.05 kcal/mol barrier gap between the hexose and pentose sulfur-addition steps collapses it,
so ~3× is carried by an ESTIMATED, unconstrained barrier rather than by mechanism), matrix
directionality — plus an experiment-prioritization machinery whose gap map is now honest.

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
| `no_verifiable_source` references | 43 | **84** (62 still feeding runtime numbers) — revised to **102 / 80 / 62 runtime** later the same day when `data/qm/` was tracked for the first time; see Open items |
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
| 6 | **SERIOUS** — Hofmann & Schieberle 1998 reports **mol %**, not ppb. The 342 / 200 ppb values require a conversion documented nowhere in this repository, on the panel's tightest contract (1.45×). | Values unchanged; a `content_verification_note` now says the derivation is unverified |

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
| Hold-out median fold error | 32.79× | 15.31× → **42.62×** after Wave O | Li 2026 wrong-row correction, then the Wave O observability refit |
| Hold-out genuine extrapolations | 0/5 pre-widening | **1/5** pre-widening (unmoved by Wave O) | Li 2026 wrong-row correction |
| Hold-out coverage at the shipped sigma | 5/8 | **4/8** | Wave O observability refit |
| Literature-row projection objective | 0.74 dex | 0.96 dex → **0.88 dex** after Wave O | Both corrections; then the Wave O refit (the 0.08 dex gain is entirely two fit-recovery rows and is **not** evidence) |
| Matrix leave-lane-out ln-sigma | 2.86 RMS on 6 residuals | **3.25** on 5 | Pratap-Singh content correction |
| Matrix leave-lane-out ln-sigma, after the Wave O refit | 3.25 on 5 | **3.25 on 5, bit-identical** | The derivation never reads the observability factors — by design |
| Test suite | 1245 passed, 0 failed | **1242 passed, 1 skipped, 2 xfailed, 0 failed** | 18 stale pins re-derived, none relaxed |

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

## What this repository is now for

A research prototype of the *architecture* alternative-protein flavor science needs —
reaction enumeration with provenance tiers, matrix-aware observability, uncertainty
propagation, and a value-of-information loop that converts every model failure into a
ranked wet-lab request — together with a documented map of the kinetic measurements the
field does not have. It is not a validated quantitative predictor, and its reports say so
on every surface.

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

| Compound | Before | After | Measured |
| --- | --- | --- | --- |
| 2-methyl-3-furanthiol (MFT) | 61.25 ppb — **5.58× under** | 235.32 ppb — **1.45× under** | 342 ppb |
| 2-furfurylthiol (FFT) | 61.44 ppb — **3.26× under** | 219.96 ppb — **1.10× over** | 200 ppb |

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
repo-wide count from 84 to 102 (see above). Two signatures, recorded in the files themselves,
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
re-anchoring the sulfur branch on new measurements, the **102 records marked
`no_verifiable_source`** whose numbers are now honestly unanchored — measured directly
across every tracked JSON and YAML file including nested records, of which **80 carry
numeric payloads** and **62 of those are consumed at runtime** — and the wet-lab campaign
([PPI/SPI protocol](docs/protocols/PPI_SPI_PRIMARY_BENCHMARK_PROTOCOL.md)) that would give
this model its first real calibration surface.
