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
   accepted norfuraneol mechanism; a fictitious hydrogen economy manufactured H₂ from
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
  was rebuilt from a verified 1994 source — which the model honestly fails at ~745×. A CI
  citation gate ([scripts/ci/citation_gate.py](scripts/ci/citation_gate.py)) now blocks
  contamination regrowth, at zero waivers.
- **Chemistry**: the real norfuraneol MFT route was implemented (the pentose≫hexose
  selectivity is now structural, not fitted); the fabricated radical chain was removed
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
| Panel | 16 benchmarks, 41 matched rows (15 literature, 26 internal-synthetic) |
| Literature rows inside 90% CI | **2/11 evaluable** (median CI width 3.0 dex) |
| Internal-synthetic rows (reproducibility only) | 18/18 — carries zero validation weight |
| External hold-out | 5/8 nominal = 3/3 in-panel re-scoring + **2/5 genuine extrapolation**; median fold error 33× |
| Strict-ready benchmarks | **0/16** |
| Reaction-chemistry lanes with generative templates | 5/16 (derived by enumeration, test-pinned) |
| Test suite | 1178 passed; 1 environmental failure (dvipng) |
| Citation gate | 0 dead DOIs, 0 waivers |

The coverage numbers *fell* during remediation — twice — and each fall is documented at
the point it happened (see the README's calibration section and the ledger). That is the
audit working as intended: each drop is the measured size of a fabricated support that was
removed. What survives at quantitative fidelity is narrow; what survives robustly is
**ordering** — pentose ≫ hexose sulfur chemistry (15.8× against a 3× floor), matrix
directionality — plus an experiment-prioritization machinery whose gap map is now honest.

## What this repository is now for

A research prototype of the *architecture* alternative-protein flavor science needs —
reaction enumeration with provenance tiers, matrix-aware observability, uncertainty
propagation, and a value-of-information loop that converts every model failure into a
ranked wet-lab request — together with a documented map of the kinetic measurements the
field does not have. It is not a validated quantitative predictor, and its reports say so
on every surface.

## Open items

The ledger's `[P]` entries are the live list. The largest: the volatile-budget
*allocation* layer (where the Wave-H refit localized the remaining sulfur residual — the
legacy heuristics that could absorb it are deliberately left dead rather than refitted),
re-anchoring the sulfur branch on new measurements, the 43 references with no verifiable
source whose numbers are now honestly unanchored, and the wet-lab campaign
([PPI/SPI protocol](docs/protocols/PPI_SPI_PRIMARY_BENCHMARK_PROTOCOL.md)) that would give
this model its first real calibration surface.
