> # ⛔ RETRACTED — 2026-08-27 (Wave S2c)
>
> **This record is retracted. Do not use it as a warrant for any shipped constant.**
>
> Its only fit target was `cys_ribose_140C_Hofmann1998`, and that benchmark's two comparator
> values are **not measurements**. Wave S2b established (~90% confidence, arithmetic exact) that
> MFT **342 ppb** and FFT **200 ppb** were derived **inside this repository**, from
> `data/benchmarks/maillard_validation_benchmarks.md` §1.3 — an abstract-reconstructed range
> table committed in `c7efbbc`, the *same commit* that created the benchmark JSON. That table's
> row reads `| Ribose + Cys, pH 5 aqueous | 140 | 30 min | ~0.02–0.05 | ~0.01–0.03 |` (MFT / FFT
> mol %). On the benchmark's declared — and itself unattested — 10 mM basis with MW 114.17:
> `0.0300 mol % × 0.010 M × 114.17 = 342.5 → 342 ppb`, and the **geometric mean** of the FFT
> band, `0.017321 mol % → 197.8 → 200 ppb`. Both targets are interior points of two invented,
> **overlapping** bands (MFT 228–571 ppb, FFT 114–342 ppb). Nothing in Hofmann & Schieberle 1998,
> nothing in Mottram & Nobrega 2002, and nothing in any retrievable literature produces 342 or 200.
>
> **Action taken.** The single constant this record moved — `thiol_addition_pentodiulose`
> 28.60 → 26.35 kcal/mol — is **reverted to 28.60** (the un-fitted `thiol_addition` class value
> Wave N shipped) in `src/barrier_constants.py`. Every other knob this record profiled was
> already reported "incumbent kept", so no other constant is affected. The benchmark's own
> `validation_contract` (1.45× / 0.09 dex) is retired in the same wave, and its `metadata.tier`
> is demoted `PRIMARY → REFERENCE`.
>
> **The cost, reported rather than buried.** On the retired benchmark itself the revert moves
> MFT **154.85 → 78.09 ppb** and FFT **267.50 → 293.67 ppb** against the fabricated 342 / 200,
> i.e. `max_ratio` **2.2086 → 4.3797** and MALE **0.2352 → 0.4041 dex**. The row got much worse.
> It has to: 26.35 existed precisely to pull MFT toward 342, and 342 is not a measurement. A
> model that looks worse against a number nobody measured has not got worse.
>
> **Same treatment as** `results/validation/hydrolysate_observability_rederivation.md`, retracted
> by Wave I when its only two fit targets turned out to be fabricated and its one applied value
> (the Methional `base_factor`) was reverted.
>
> **What this does *not* retract.** The fit-target accounting is unchanged.
> `cys_ribose_140C_Hofmann1998` is still a declared fit target of two other **live** records —
> `projection_constant_refit.json` and `sulfur_barrier_refit_hofmann.json` — so it stays flagged
> `fitted_row: true` in `prediction_uncertainty.json` and stays **out** of the honest
> literature-coverage numerator and denominator. That is the correct outcome: a value with no
> verifiable source must not count as literature evidence in *either* direction.
> `scripts/ci/fit_target_gate.py` skips retracted records, so the requirement on this benchmark
> now flows from those two records rather than three; the flag it enforces is unchanged and the
> gate is not weakened.
>
> **Do not re-run** `scripts/generators/refit_thiol_addition_pentodiulose_hofmann.py` against
> this benchmark. It needs a real target first — see the ILL pack in
> `tasks/audit_remediation.md` "## Wave S2b" §(f), then a rebuild in native **mol %**.
>
> **THE SULFUR BRANCH HAS ZERO ABSOLUTE LITERATURE ANCHORS.**
>
> The text below is preserved **verbatim as the forensic record of what was done**, not as a
> current claim. In particular its first line, calling this benchmark "the ONLY surviving
> literature constraint on the sulfur branch", is **false** — it was never a literature
> constraint at all.

---

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
