# Projection constant refit (v2 Arrhenius budget)

Objective: `mean |log10(predicted_ppb / measured_ppb)| over literature-sourced benchmark rows`

Fit targets: 12 literature-sourced benchmarks (24/24 rows scored).

Synthetic (Internal2026 / ProtocolPilot2026) rows and the `external_validation/` hold-out are excluded by assertion, not convention.

## Shipped vs fitted (2026-08-27, Wave H)

**This fit is NOT applied to the runtime.** The script has never written to `src/projection.py`; after the Wave G1 chemistry rebuild the shipped value and the refit optimum are no longer the same number, so both are stated here.

| Constant | Shipped (runtime) | Refit optimum |
| --- | --- | --- |
| baseline_volatile_yield_fraction | 1.000e-06 | 1.000e-06 |
| reference_conversion_time_min | 1.259e+04 | 2512 |
| objective (all literature rows) | 0.8777 dex | 0.7543 dex |

**Decision: FIT NOT APPLIED — incumbent reference_conversion_time_min retained.**

tau_ref is a SINGLE GLOBAL SCALE on the volatile budget. Moving it does not improve any mechanism; it trades one lane against another, and this panel's mean-|log10| objective takes that trade happily. At the refit optimum the Hofmann1998 MFT residual collapses to 0.048 dex (from 0.75) while the Resconi 2023 furfural residual grows to 1.281 dex, i.e. ~19x OVER. The arithmetic is decisive: at the shipped value the total volatile budget at the Hofmann conditions is ~1050 ppb against a measured MFT+FFT of 542 ppb, so the budget is already the right order of magnitude and the sulfur deficit is an ALLOCATION deficit (furfural, unmeasured in that benchmark, takes ~78% of the pool) — see results/validation/sulfur_barrier_refit_hofmann.md. Raising the global budget 5x to close a 5.6x allocation gap means supplying ~5250 ppb of volatiles to explain 542 ppb of measured product and sending the other ~4700 ppb to the species nobody measured in that experiment. That is not a calibration, it is a way of making the Wave H finding invisible. The incumbent is therefore retained and the gap is reported.

What would justify applying it: A refit of the ALLOCATION (the depth-bias and direct-sulfur heuristics, both unconstrained legacy fits from the quarantined-target era), done first and with its own literature constraints, after which a global scale would be fitting a scale rather than absorbing a selectivity error.

## Fitted constants

| Constant | Value |
| --- | --- |
| baseline_volatile_yield_fraction | 1.000e-06 |
| reference_conversion_time_min | 2512 |

## Fixed (not fitted)

| Constant | Value | Basis |
| --- | --- | --- |
| apparent_activation_energy_kj_mol | 120.0 | data/lit/arrhenius_params.yml `enolisation` |
| reference_temperature_kelvin | 423.15 | 150 C panel anchor |
| conversion_ceiling_fraction | 1.0 | mass conservation |

Objective at fit: **0.7543 dex** over all 24 literature rows; **0.8074 dex** over the 11 rows the projection layer actually moves (the rest run on the matrix-only and safety lanes, which never touch the volatile budget).

## Versus v1 (severity sigmoid + double Boltzmann)

| Objective | v1 | v2 (this fit) |
| --- | --- | --- |
| All literature rows | 0.1469 | 0.7543 |
| Projection-sensitive rows | 0.3089 | 0.8074 |

The v1 objective is BETTER. Essentially all of the gap is the de-duplication of the Boltzmann factor, not the budget retune: the budget retune is a single global scale that the tau_ref fit absorbs, whereas removing the unphysical T/2.19 sharpening changes how the budget is split between competing species. In other words, ~0.27 dex of the panel's previous agreement on multi-analyte free-precursor benchmarks was being carried by a selectivity term evaluated at less than half the physical temperature. That is a finding, not a regression to be tuned away.

## Per-row residuals at the fit

| Benchmark | Compound | Measured ppb | Predicted ppb | \|log10 ratio\| | Projection-sensitive |
| --- | --- | --- | --- | --- | --- |
| cml_cel_commercial_pbma_Foods2023 | Nε-(Carboxymethyl)lysine (CML) | 3.2e+04 | 26.59 | 3.081 | no |
| cml_cel_commercial_pbma_Foods2023 | Nε-(Carboxyethyl)lysine (CEL) | 5.5e+04 | 51.68 | 3.027 | no |
| furosine_extrusion_crossover_140C_RamirezJimenez2000 | furosine | 1.74e+04 | 86.55 | 2.303 | no |
| thiamine_cys_glucose_120C_Bolton1994 | 2-Methyl-3-furanthiol (MFT) | 13 | 0.08721 | 2.173 | yes |
| resconi_2023_pbma_beef_identity_benchmark | furfural | 1040 | 1.988e+04 | 1.281 | yes |
| spi_hvp_xylose_120C_PMC9905368 | 2-Furfurylthiol (FFT) | 0.42 | 0.02214 | 1.278 | yes |
| wheat_gluten_hvp_xylose_120C_PMC9905368 | 2-Furfurylthiol (FFT) | 0.61 | 0.06502 | 0.972 | yes |
| spi_hvp_xylose_120C_PMC9905368 | 2-Methyl-3-furanthiol (MFT) | 0.18 | 0.02679 | 0.827 | yes |
| acrylamide_spi_extrusion_130C_ACSRef3 | acrylamide | 62.62 | 9.748 | 0.808 | no |
| wheat_gluten_hvp_xylose_120C_PMC9905368 | Methional | 3.44 | 17.26 | 0.700 | yes |
| wheat_gluten_hvp_xylose_120C_PMC9905368 | 2-Methyl-3-furanthiol (MFT) | 0.34 | 0.07869 | 0.636 | yes |
| spi_hvp_xylose_120C_PMC9905368 | Methional | 1.88 | 7.214 | 0.584 | yes |
| thiamine_cys_xylose_145C_Cerny2008 | 2-Methyl-3-furanthiol (MFT) | 2.47 | 3.863 | 0.194 | yes |
| cys_ribose_140C_Hofmann1998 | 2-furfurylthiol | 200 | 307 | 0.186 | yes |
| cys_ribose_140C_Hofmann1998 | 2-methyl-3-furanthiol | 342 | 306.1 | 0.048 | yes |
| pea_isolate_40C_PratapSingh2021 | hexanal | 260 | 260.6 | 0.001 | no |
| pea_isolate_40C_PratapSingh2021 | hexanol | 80 | 80.1 | 0.001 | no |
| soy_isolate_40C_PratapSingh2021 | hexanol | 120 | 119.9 | 0.000 | no |
| pea_isolate_40C_PratapSingh2021 | 2-pentylfuran | 638 | 638.3 | 0.000 | no |
| soy_isolate_40C_PratapSingh2021 | hexanal | 380 | 379.9 | 0.000 | no |
| soy_isolate_40C_PratapSingh2021 | 2-pentylfuran | 2492 | 2492 | 0.000 | no |
| pea_isolate_uht_140C_Trikusuma2019 | hexanal | 782 | 782 | 0.000 | no |
| pea_isolate_uht_140C_Trikusuma2019 | 2-pentylfuran | 163 | 163 | 0.000 | no |
| pea_isolate_uht_140C_Trikusuma2019 | nonanal | 24 | 24 | 0.000 | no |
