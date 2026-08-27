# Sulfur-branch barrier refit against Hofmann1998 (Wave H)

Fit target: `cys_ribose_140C_Hofmann1998.json` (DOI 10.1021/jf9705983) — the ONLY surviving literature constraint on the sulfur branch after the Mottram1994 / Farmer1999 quarantine.

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
| `thiol_addition_norfuraneol` | [23.30, 29.65] | H2S/thiol addition to a carbonyl; low = thiohemiacetal_formation (23.30), high = thiol_addition_hexose (29.65), both already in FAST_BARRIERS | 26.85 | 0.0789 | IMMATERIAL — best achievable gain 0.0099 dex < 0.02 dex; incumbent kept |
| `furanone_cyclisation` | [26.00, 30.00] | cyclodehydration class; dehydration / 2,3-enolisation / furanone_formation are all 28.0 in FAST_BARRIERS, taken +/- 2 kcal/mol | 28.00 | 0.0690 | IMMATERIAL — best achievable gain 0.0000 dex < 0.02 dex; incumbent kept |
| `thiohemiacetal_formation` | [21.00, 26.80] | thiohemiacetal formation must stay below its own dehydration step, thiol_dehydration (26.80) | 23.30 | 0.0000 | NOT IDENTIFIABLE — flat profile; incumbent kept |
| `thiol_dehydration` | [23.30, 28.00] | bounded below by thiohemiacetal_formation (23.30) and above by the generic dehydration barrier (28.00) | 26.80 | 0.0407 | IMMATERIAL — best achievable gain 0.0096 dex < 0.02 dex; incumbent kept |
| `thiol_addition` | reported, not fitted | After Wave G1 this key labels only `Thiol_Addition_Legacy_Shortcut` (the demoted fabricated one-step MFT route) and the lumped `Thiol_Addition_H2` step. Fitting it would tune the route this wave exists to retire; its derivative is reported instead. | 28.60 | 0.0047 | NOT FITTED |

Objective: **0.6298 dex** (incumbent) -> **0.6298 dex** (adopted).

| Compound | Measured ppb | Predicted ppb (incumbent) | Predicted ppb (adopted) | Fold error (adopted) |
| --- | --- | --- | --- | --- |
| 2-methyl-3-furanthiol | 342 | 61.25 | 61.25 | 5.58x under |
| 2-furfurylthiol | 200 | 61.44 | 61.44 | 3.26x under |

## What this fit cannot fix

The profiles saturate well inside the defensible ranges: the MFT ceiling is set by the shared upstream span (Schiff/Amadori) and by the volatile-budget allocation, in which furfural — unmeasured in this benchmark — takes ~78% of a total budget that is itself the right order of magnitude (~1050 ppb against a measured MFT+FFT of 542 ppb). No barrier value in any defensible range reproduces the measured absolute sulfur yields; the remaining factor of ~5 on MFT and ~3 on FFT is an allocation deficit, not a barrier deficit. That is the finding this refit exists to produce, and it is reported rather than tuned away.

With every fitted knob pinned at the bottom of its defensible range the objective is still **0.6102 dex**.
