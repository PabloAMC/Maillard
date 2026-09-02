> # ⛔ RETRACTED — 2026-08-27 (Wave I)
>
> **This record is retracted. Do not use it as a warrant for any shipped constant.**
>
> Its only fit targets were `spi_hvp_xylose_120C_PMC9905368` and
> `wheat_gluten_hvp_xylose_120C_PMC9905368`. The cold-start red team established that both are
> **fabricated**: the cited paper `10.1007/s10068-022-01194-w` reacts protein hydrolysates with
> **glucose and fructose at pH 7.5 for 90 min** and reports only **relative GC-MS peak areas**.
> It never mentions 2-furfurylthiol or 2-methyl-3-furanthiol and contains no absolute
> concentration for any analyte, so the `conc_ppb` values fitted against here have **no possible
> source**. Both files are now in `data/benchmarks/quarantined/`.
>
> **Action taken.** The single constant this record moved — Methional `base_factor`
> 0.0045 → 0.05623 — is **reverted** in
> `src/recommend.py::_HYDROLYSATE_SULFUR_OBSERVABILITY_PROFILES`. The other three compounds were
> not moved here (saturated / not derivable), so no other constant is affected.
>
> **Why it matters beyond the constant.** The two Methional rows this fit produced were then
> *scored* as literature coverage in the headline — a value fitted to a benchmark and then
> reported as agreement with that same benchmark. `scripts/ci/fit_target_gate.py` (Wave I) now
> makes that circularity a blocking CI failure.
>
> The text below is preserved **verbatim as the forensic record of what was done**, not as a
> current claim.

---

# Hydrolysate sulfur observability factors — re-derivation (Wave H)

`src/recommend.py::_HYDROLYSATE_SULFUR_OBSERVABILITY_PROFILES` holds one `base_factor` per compound: the fraction of the modelled volatile taken to survive to the headspace in an HVP matrix. These constrain the PRODUCT of the lane, not a barrier, so after Wave G1's chemistry change and Wave H's Hofmann-only barrier refit they are re-derived here — the same way they were originally derived, against the same two literature benchmarks — rather than being compensated for by bending a barrier.

Fit targets: spi_hvp_xylose_120C_PMC9905368.json, wheat_gluten_hvp_xylose_120C_PMC9905368.json (literature only; the synthetic lanes and the hold-out are excluded by assertion).

Constraint: base_factor is a surviving fraction and lives in [0.0001, 1.0] — the same clamp `_resolve_upstream_observability_factor` applies at runtime.

## What this re-derivation changed (2026-08-27)

| Compound | From | To | Objective | Residuals |
| --- | --- | --- | --- | --- |
| Methional | 0.0045 | 0.05623 | 1.1560 -> 0.0619 dex | 16.2x under (SPI) / 12.4x under (wheat gluten) -> 1.30x under (SPI) / 1.01x (wheat gluten) |

This script is idempotent: re-run after the value is applied and it reports no further move, which is the convergence check. The table below is that re-run, so its residuals are the CURRENT ones.

| Compound | Incumbent | Unconstrained optimum | Decision |
| --- | --- | --- | --- |
| Methional | 0.05623 | 0.06391 | IMMATERIAL — best achievable gain 0.0000 dex < 0.02 dex; incumbent kept |
| 2-Furfurylthiol | 0.13 | 8.65 | SATURATED — the unconstrained optimum is 8.65, 8.6x ABOVE the physical maximum of 1.0 for a surviving fraction. This layer cannot explain the residual (it can only suppress, and the model already under-predicts this lane by 66.5x). Incumbent kept. |
| 2-Methyl-3-furanthiol | 0.13 | 3.493 | SATURATED — the unconstrained optimum is 3.49, 3.5x ABOVE the physical maximum of 1.0 for a surviving fraction. This layer cannot explain the residual (it can only suppress, and the model already under-predicts this lane by 26.9x). Incumbent kept. |
| bis(2-methyl-3-furyl) disulfide | 0.18 | n/a | NOT DERIVABLE — this compound has no row in any LITERATURE benchmark of this lane; its only comparators are the synthetic Internal2026 / ProtocolPilot2026 snapshots, which are forbidden as fit targets. The incumbent is kept as an unconstrained legacy estimate. |

## Per-row residuals at the incumbent

| Benchmark | Compound | Measured ppb | Predicted ppb | Fold error |
| --- | --- | --- | --- | --- |
| spi_hvp_xylose_120C_PMC9905368 | Methional | 1.88 | 1.446 | 1.3x under |
| wheat_gluten_hvp_xylose_120C_PMC9905368 | Methional | 3.44 | 3.461 | 1.0x over |
| spi_hvp_xylose_120C_PMC9905368 | 2-Furfurylthiol (FFT) | 0.42 | 0.004439 | 94.6x under |
| wheat_gluten_hvp_xylose_120C_PMC9905368 | 2-Furfurylthiol (FFT) | 0.61 | 0.01304 | 46.8x under |
| spi_hvp_xylose_120C_PMC9905368 | 2-Methyl-3-furanthiol (MFT) | 0.18 | 0.005372 | 33.5x under |
| wheat_gluten_hvp_xylose_120C_PMC9905368 | 2-Methyl-3-furanthiol (MFT) | 0.34 | 0.01578 | 21.5x under |
