# Sulfur-branch barrier refit against Hofmann1998 (Wave H)

Fit target: `cys_ribose_140C_Hofmann1998.json` (DOI 10.1021/jf9705983) — the ONLY surviving literature constraint on the sulfur branch after the Mottram1994 / Farmer1999 quarantine.

> **CORRECTION — 2026-08-27 (Wave S2c). The sentence directly above is false in the word
> "literature", and is kept verbatim because it is what the repo believed.** Wave S2b traced this
> benchmark's MFT **342 ppb** / FFT **200 ppb** to `data/benchmarks/maillard_validation_benchmarks.md`
> §1.3 — an abstract-reconstructed range table committed in `c7efbbc`, the *same commit* that
> created the benchmark JSON — whose row reads `~0.02–0.05` mol % (MFT) and `~0.01–0.03` mol %
> (FFT). On the file's declared 10 mM basis with MW 114.17: `0.0300 mol % → 342.5 → 342 ppb`, and
> the geometric mean of the FFT band `0.017321 mol % → 197.8 → 200 ppb`. Both are interior points
> of two invented, overlapping bands (~90% confidence, arithmetic exact). **This record's fit
> target was never a literature constraint. THE SULFUR BRANCH HAS ZERO ABSOLUTE LITERATURE
> ANCHORS.**
>
> **Not retracted, unlike its companion.** The constant this record moved,
> `thiol_addition_norfuraneol` = 26.85, sits on a family **no step emits** since Wave N removed
> the norfuraneol → MFT step on isotope evidence, so reverting it would move no prediction while
> making the provenance record harder to read. It is annotated instead. The companion record
> `results/validation/sulfur_barrier_refit_pentodiulose.{json,md}` **is** retracted and its
> constant reverted (`thiol_addition_pentodiulose` 26.35 → 28.60), because that one is live on the
> shipping network. Keeping *this* record live also keeps `cys_ribose_140C_Hofmann1998` a declared
> fit target for `scripts/ci/fit_target_gate.py`, hence `fitted_row: true` and hence **out** of the
> honest literature-coverage numerator and denominator — the right place for a value with no
> verifiable source. The gate is not weakened.

Objective: `mean |log10(predicted_ppb / measured_ppb)| over the matched rows of cys_ribose_140C_Hofmann1998`

Synthetic (Internal2026 / ProtocolPilot2026) lanes and the `external_validation/` hold-out are excluded **by assertion**, not convention. The xylose HVP lanes are excluded because they are constrained by their own per-lane `upstream_observability_factor`, re-derived separately.

## Decision rules

1. **Identifiability** — a profile flat to 1e-06 dex over the whole defensible range keeps its incumbent.
2. **Materiality** — a knob moves only for a gain of at least 0.02 dex.
3. **Conservative edge** — among values within 0.01 dex of the profile minimum, the LARGEST barrier is adopted.

## What this refit changed (2026-08-27)

| Knob | From | To | Objective | Hofmann MFT | Hofmann FFT |
| --- | --- | --- | --- | --- | --- |
| `thiol_addition_norfuraneol` | 28.60 | 26.85 | 0.6987 -> 0.6298 dex | 7.83x -> 5.58x under | 3.19x -> 3.26x under |

This script is idempotent: re-run after the fit is applied and it reports no further move, which is the convergence check. The section below is that re-run.

## Result

| Knob | Range | Basis | Incumbent | Profile span (dex) | Decision |
| --- | --- | --- | --- | --- | --- |
| `thiol_addition_norfuraneol` | [23.30, 29.65] | H2S/thiol addition to a carbonyl; low = thiohemiacetal_formation (23.30), high = thiol_addition_hexose (29.65), both already in FAST_BARRIERS | 26.85 | 0.3653 | IMMATERIAL — best achievable gain 0.0183 dex < 0.02 dex; incumbent kept |
| `furanone_cyclisation` | [26.00, 30.00] | cyclodehydration class; dehydration / 2,3-enolisation / furanone_formation are all 28.0 in FAST_BARRIERS, taken +/- 2 kcal/mol | 28.00 | 0.5152 | MOVED 28.00 -> 26.65 (conservative edge of the 26.00-26.65 indifference band; gain 0.0730 dex) |
| `thiohemiacetal_formation` | [21.00, 26.80] | thiohemiacetal formation must stay below its own dehydration step, thiol_dehydration (26.80) | 23.30 | 0.0201 | IMMATERIAL — best achievable gain 0.0199 dex < 0.02 dex; incumbent kept |
| `thiol_dehydration` | [23.30, 28.00] | bounded below by thiohemiacetal_formation (23.30) and above by the generic dehydration barrier (28.00) | 26.80 | 0.0450 | MOVED 26.80 -> 27.70 (conservative edge of the 27.25-27.70 indifference band; gain 0.0242 dex) |
| `thiol_addition` | reported, not fitted | After Wave G1 this key labels only `Thiol_Addition_Legacy_Shortcut` (the demoted fabricated one-step MFT route) and the lumped `Thiol_Addition_H2` step. Fitting it would tune the route this wave exists to retire; its derivative is reported instead. | 28.60 | 0.1011 | NOT FITTED |

Objective: **0.1018 dex** (incumbent) -> **0.0558 dex** (adopted).

| Compound | Measured ppb | Predicted ppb (incumbent) | Predicted ppb (adopted) | Fold error (adopted) |
| --- | --- | --- | --- | --- |
| 2-methyl-3-furanthiol | 342 | 235.3 | 301.7 | 1.13x under |
| 2-furfurylthiol | 200 | 220 | 175.4 | 1.14x under |

## What this fit cannot fix

The profiles saturate well inside the defensible ranges: the MFT ceiling is set by the shared upstream span (Schiff/Amadori) and by the volatile-budget allocation, in which furfural — unmeasured in this benchmark — takes ~78% of a total budget that is itself the right order of magnitude (~1050 ppb against a measured MFT+FFT of 542 ppb). No barrier value in any defensible range reproduces the measured absolute sulfur yields; the remaining factor of ~5 on MFT and ~3 on FFT is an allocation deficit, not a barrier deficit. That is the finding this refit exists to produce, and it is reported rather than tuned away.

With every fitted knob pinned at the bottom of its defensible range the objective is still **0.0282 dex**.


## NOT APPLIED

**The adopted point above is NOT what ships.** This run was report-only. Shipped vs adopted:

| Constant | Shipped | Adopted here |
| --- | ---: | ---: |
| `furanone_cyclisation` | 28.00 | 26.65 |
| `thiol_dehydration` | 26.80 | 27.70 |

Wave I (2026-08-27) re-ran this record because the network changed underneath it: decoupling the MFT step from the pyrazine-supplied hydrogen pool took cys_ribose_140C_Hofmann1998 from 5.58x under to 1.45x under with NO barrier changed. Applying a barrier refit on top of a chemistry change, in the same pass, would entangle the two and make it impossible to say afterwards which one produced the agreement. The record is written; the constants are not moved. Applying it is an OPEN OWNER ITEM.