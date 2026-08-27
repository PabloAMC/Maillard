# `thiol_addition_pentodiulose` refit against Hofmann1998 (Wave P item 1)

Fit target: `cys_ribose_140C_Hofmann1998.json` (DOI 10.1021/jf9705983) — the ONLY surviving literature constraint on the sulfur branch.

## The caveat this number carries

The fit target's own `content_verification_note` says, verbatim:

> The paper's abstract reports MFT/FFT yields in mol % (e.g. MFT 1.4 mol % for its best precursor system), not ppb. The 342/200 ppb values here require a mol%->ppb conversion (system volume, precursor moles, molecular weights) that is NOT documented anywhere in this repo, and the full text is paywalled so the values could not be confirmed against the paper's tables. This is the panel's tightest contract (1.45x / 0.09 dex) resting on an unverified derivation. Needs: full-text access to 10.1021/jf9705983 and a written conversion, or replacement with the paper's native mol % as the target unit.

This constant is fitted against an anchor whose mol%->ppb conversion is UNVERIFIED. If the conversion is wrong, the resulting error is LOCALISED in this one barrier rather than distributed across the route. That is the argument for doing the fit at all, and it is also the reason the value must never be quoted as a measured or literature-derived barrier.

## Result

| Knob | Range | Incumbent | Profile span (dex) | argmin | Decision |
| --- | --- | ---: | ---: | ---: | --- |
| `thiol_addition_pentodiulose` | [23.30, 29.65] | 26.35 | 0.3224 | 23.30 | IMMATERIAL — best achievable gain 0.0099 dex < 0.02 dex; incumbent kept |

Objective: **0.0935 dex** (incumbent) -> **0.0935 dex** (adopted).

Boundary check: `hit_a_bound = True`. The objective minimum sits AT a range endpoint. The data want a value the table's own sulfur-addition class envelope does not contain, i.e. the profile SATURATES: the residual is not removable by this barrier. The conservative-edge rule is what keeps the adopted value off the boundary.

| Compound | Measured ppb | Predicted (incumbent) | Predicted (adopted) | Fold error (adopted) |
| --- | ---: | ---: | ---: | --- |
| 2-methyl-3-furanthiol | 342 | 242.4 | 242.4 | 1.4110x under |
| 2-furfurylthiol | 200 | 218 | 218 | 1.0900x over |

## FFT co-movement (measured, not assumed)

The fitted knob sits on the MFT lane only. FFT shares the upstream sugar flux, so it can move without being fitted. Measured over the same grid.

Over the whole search range MFT spans **3.534x** and FFT spans **1.249x**. FFT moves OPPOSITE to MFT: lowering the MFT barrier diverts shared upstream flux into the MFT lane and FFT falls. The two rows are therefore NOT independent evidence, which is exactly why one knob is fitted and not two.

## Reported, not fitted

| Knob | Incumbent | Objective span over +/-2 kcal (dex) | Why not fitted |
| --- | ---: | ---: | --- |
| `furanone_cyclisation` | 28.00 | 0.0000 | The competing branch of the same 1-deoxyosone fork (-> norfuraneol). Wave N set it EQUAL to `deoxyosone_reduction` (both 28.0) precisely to express no prior preference at that fork; fitting it here would encode a preference derived from two rows of one benchmark. |
| `deoxyosone_reduction` | 28.00 | 0.3850 | The upstream reduction that feeds the fitted step. It is in series with it, so on a single benchmark the two are not separately identifiable -- fitting both would split one measurable quantity across two constants and report the split as knowledge. |
| `thiohemiacetal_formation` | 23.30 | 0.0033 | FFT-side, not MFT-side. It moves the OTHER measured row of this benchmark, so fitting it would let the fit trade FFT accuracy for MFT accuracy across the two rows -- the classic two-row/two-knob recovery. |
| `thiol_dehydration` | 26.80 | 0.0285 | Same reason as `thiohemiacetal_formation`: FFT-side. |
| `thiol_addition` | 28.60 | 0.1011 | After Wave G1/Wave I this key labels only the DEMOTED hexose legacy shortcut and the lumped `Thiol_Addition_H2` FFT step. Fitting it would tune a lump the campaign exists to retire. |

## What this fit cannot fix

With the knob pinned at the bottom of its defensible range (23.30 kcal/mol) the objective is still **0.0836 dex**.

## Fit leverage

ONE free barrier against 2 rows is 0.50 parameters per row, at src.fit_target_index.FIT_LEVERAGE_THRESHOLD (0.5). cys_ribose_140C_Hofmann1998 therefore stays OUT of the honest literature-coverage numerator and denominator and is reported in the fitted-row bucket. That classification is UNCHANGED by this refit -- the benchmark was already a declared fit target of the Wave H record. What changes is which constant carries the fit, and how much of the residual it can absorb.
