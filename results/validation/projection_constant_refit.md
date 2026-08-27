# Projection constant refit (v2 Arrhenius budget)

Objective: `mean |log10(predicted_ppb / measured_ppb)| over literature-sourced benchmark rows`

Fit targets: 10 literature-sourced benchmarks (16/16 rows scored).

Synthetic (Internal2026 / ProtocolPilot2026) rows and the `external_validation/` hold-out are excluded by assertion, not convention.

## Shipped vs fitted (2026-08-27, Wave H)

**This fit is NOT applied to the runtime.** The script has never written to `src/projection.py`; after the Wave G1 chemistry rebuild the shipped value and the refit optimum are no longer the same number, so both are stated here.

| Constant | Shipped (runtime) | Refit optimum |
| --- | --- | --- |
| baseline_volatile_yield_fraction | 1.000e-06 | 1.000e-06 |
| reference_conversion_time_min | 1.259e+04 | 1e+04 |
| objective (all literature rows) | 0.8879 dex | 0.8817 dex |

**Decision: FIT NOT APPLIED — incumbent reference_conversion_time_min retained.**

tau_ref is a SINGLE GLOBAL SCALE on the volatile budget. Moving it does not improve any mechanism; it trades one lane against another, and this panel's mean-|log10| objective takes that trade happily. At the refit optimum the Hofmann1998 MFT residual collapses to 0.0497 dex (from 0.1495 dex shipped) while the Resconi 2023 furfural residual grows to 0.7678 dex (from 0.6681 dex shipped), i.e. from 4.7x to 5.9x OVER. The arithmetic is decisive: at the shipped value the total predicted pool over the five species this benchmark's system tracks is 1152.3 ppb against a measured MFT+FFT of 542 ppb (derived live, furfural 511.9, 2-methyl-3-furanthiol 242.4, 2-furfurylthiol 218.0, bis(2-methyl-3-furyl) disulfide 120.6, 2,5-dimethylpyrazine 59.4), so the budget is already the right order of magnitude and the sulfur deficit is an ALLOCATION deficit — furfural, unmeasured in that benchmark, takes 44% of the pool. See results/validation/sulfur_barrier_refit_pentodiulose.md (the Wave P refit of `thiol_addition_pentodiulose` against this same benchmark; it SUPERSEDES sulfur_barrier_refit_hofmann.md, which is stale because it profiles `thiol_addition_norfuraneol`, a family the isotope evidence retired). Raising the global budget 1.26x to close the current MFT gap (1.41x under) means supplying ~1451 ppb of volatiles at these conditions to explain 542 ppb of measured product, sending the rest to species nobody measured, and pushing FFT from 1.09x over further over. That is not a calibration, it is a way of making the allocation finding invisible. The incumbent is therefore retained and the gap is reported.

What would justify applying it: A refit of the ALLOCATION (the depth-bias and direct-sulfur heuristics, both unconstrained legacy fits from the quarantined-target era), done first and with its own literature constraints, after which a global scale would be fitting a scale rather than absorbing a selectivity error.

## Fitted constants

| Constant | Value |
| --- | --- |
| baseline_volatile_yield_fraction | 1.000e-06 |
| reference_conversion_time_min | 1e+04 |

## Fixed (not fitted)

| Constant | Value | Basis |
| --- | --- | --- |
| apparent_activation_energy_kj_mol | 120.0 | data/lit/arrhenius_params.yml `enolisation` |
| reference_temperature_kelvin | 423.15 | 150 C panel anchor |
| conversion_ceiling_fraction | 1.0 | mass conservation |

Objective at fit: **0.8817 dex** over all 16 literature rows; **0.8285 dex** over the 5 rows the projection layer actually moves (the rest run on the matrix-only and safety lanes, which never touch the volatile budget).

## Versus v1 (severity sigmoid + double Boltzmann)

| Objective | v1 | v2 (this fit) |
| --- | --- | --- |
| All literature rows | 0.1469 | 0.8817 |
| Projection-sensitive rows | 0.3089 | 0.8285 |

The v1 objective is BETTER. Essentially all of the gap is the de-duplication of the Boltzmann factor, not the budget retune: the budget retune is a single global scale that the tau_ref fit absorbs, whereas removing the unphysical T/2.19 sharpening changes how the budget is split between competing species. In other words, ~0.27 dex of the panel's previous agreement on multi-analyte free-precursor benchmarks was being carried by a selectivity term evaluated at less than half the physical temperature. That is a finding, not a regression to be tuned away.

## Per-row residuals at the fit

| Benchmark | Compound | Measured ppb | Predicted ppb | \|log10 ratio\| | Projection-sensitive |
| --- | --- | --- | --- | --- | --- |
| cml_cel_commercial_pbma_Foods2023 | Nε-(Carboxymethyl)lysine (CML) | 3.2e+04 | 26.59 | 3.081 | no |
| cml_cel_commercial_pbma_Foods2023 | Nε-(Carboxyethyl)lysine (CEL) | 5.5e+04 | 51.68 | 3.027 | no |
| thiamine_cys_glucose_120C_Bolton1994 | 2-Methyl-3-furanthiol (MFT) | 13 | 0.02142 | 2.783 | yes |
| furosine_extrusion_crossover_140C_RamirezJimenez2000 | furosine | 1.74e+04 | 86.55 | 2.303 | no |
| acrylamide_spi_extrusion_130C_ACSRef3 | acrylamide | 150 | 9.748 | 1.187 | no |
| resconi_2023_pbma_beef_identity_benchmark | furfural | 715.2 | 4191 | 0.768 | yes |
| thiamine_cys_xylose_145C_Cerny2008 | 2-Methyl-3-furanthiol (MFT) | 2.47 | 0.9728 | 0.405 | yes |
| pea_isolate_uht_140C_Trikusuma2019 | nonanal | 24 | 10.56 | 0.357 | no |
| cys_ribose_140C_Hofmann1998 | 2-furfurylthiol | 200 | 274.3 | 0.137 | yes |
| cys_ribose_140C_Hofmann1998 | 2-methyl-3-furanthiol | 342 | 305 | 0.050 | yes |
| soy_isolate_40C_PratapSingh2021 | hexanal | 1622 | 1640 | 0.005 | no |
| pea_isolate_40C_PratapSingh2021 | hexanal | 1138 | 1125 | 0.005 | no |
| pea_isolate_40C_PratapSingh2021 | 2-pentylfuran | 638 | 638.3 | 0.000 | no |
| soy_isolate_40C_PratapSingh2021 | 2-pentylfuran | 2492 | 2492 | 0.000 | no |
| pea_isolate_uht_140C_Trikusuma2019 | hexanal | 782 | 782 | 0.000 | no |
| pea_isolate_uht_140C_Trikusuma2019 | 2-pentylfuran | 163 | 163 | 0.000 | no |
