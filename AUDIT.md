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
| Literature rows inside 90% CI | **1/3 evaluable** (median CI width **0.95 dex** — it widened from 0.85 under Wave S1's additive propagator *without* buying any coverage), with fitted rows removed from numerator *and* denominator |
| Fitted rows (constants back-solved from the benchmark) | 2/2 — **not evidence**; both would previously have been counted as literature hits |
| Internal-synthetic rows (reproducibility only) | 18/18 — carries zero validation weight |
| Benchmarks without blocking gaps | **0/6 predictive** (+ **0/4** fit-recovery, 4/4 synthetic; the **4/14** aggregate is all non-evidence — this row still read 1/4 and 5/14 after Wave P moved it, and was corrected in Wave S1). Fit-recovery fell 3/4 → 1/4 in Round 3 when the two Pratap-Singh benchmarks were corrected against the paper and stopped recovering. Wave O then refitted their constants onto the verified values and they recover again — at `pass-no-ranking`, which is not `pass`, so **none of these counts moved**. A refit changes the size of a recovery, never its evidential status. |
| External hold-out | **1/5 genuine extrapolations** at the pre-widening prior — and 1/5 under the wider one too since Wave O; median fold error **93.68×**, worst **2474×**, coverage **3/8** at the shipped sigma. Only 4 of the 8 points are measurements at all. The median was 32.79× until Round 3 corrected two wrong table rows in the Li 2026 bundle (→ 15.31×, predictions unmoved), rose to 42.62× when Wave O refitted the ambient hexanal observability onto the paper's verified anchor, and rose again to 93.68× when Wave R found the Liu hold-out's two reference values matched nothing in their source and replaced them with the thesis's own Table 2.7. **Every rise came from correcting data toward the literature, and predictions moved in none of them**; see Round 3, Wave O and Wave R. |
| Strict-ready benchmarks | **0/14** |
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
   it was removed — the accepted mechanism names the reductant and it is the amino acid
   (Blank & Fay 1996, `10.1021/jf950439o`; Kerler et al. 2010, `10.1002/9781444317770.ch3`).
   Measured: disabling aminoketone condensation now leaves the DMHF step standing (it used to
   take it to zero).

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
