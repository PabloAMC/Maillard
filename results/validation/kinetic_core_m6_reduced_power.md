# Case M6: falsification arithmetic and reduced-model power check (2026-11-07)

Case: `docs/reference/cases/M6_pyrazines.md (branch m6-case-ready)`. Model: `src/m6_reduced_power.py`. SYNTHETIC data for job (c); job (a) is arithmetic on the spent / seen printed values. No engine run or fit; no decisive, reserved or TEST value read. Judge: the adopted sec. 7 item 5 / M1.1 (`src/kinetic_core/predictive_densities.py`).

## Job (a): falsifications on data on disk (falsify-only; spent data never credit)

**FX-C1, Zhou 2024's six zero-order R2 vs a free-aminoketone pool with Elfsson's class and no first-order exit.** p = max over noise SD and form of P(R2 >= printed):

| unit | tau (h) | R2 noiseless | p max | p Holm | k2 needed for a 10-min lag / Elfsson's fastest corner (dex) |
|---|---:|---:|---:|---:|---:|
| PZ_100C | 21.6 | 0.821 | 0.0115 | 0.0231 | 4.22 |
| PZ_110C | 10.5 | 0.821 | 0.0027 | 0.0135 | 3.60 |
| PZ_120C | 6.3 | 0.822 | 0.0260 | 0.0260 | 3.15 |
| DMP_100C | 60.9 | 0.821 | 0.0043 | 0.0172 | 5.13 |
| DMP_110C | 29.4 | 0.821 | 0.0015 | 0.0090 | 4.49 |
| DMP_120C | 16.0 | 0.821 | 0.0069 | 0.0205 | 3.97 |

Falsified: **True**. Smallest first-order AK exit that survives (Holm, every unit): **ks >= 0.0281 /min** at 20 mM MGO. At that edge (100 C, pH 8, Elfsson's fastest corner): AK 51 uM, Strecker flux 1.44 uM/min, kS 3.59e-03 M^-1 min^-1 = B18 + 2.08 dex = Hofmann-implied -0.22 dex; DMP per Strecker event 0.0049; exit product 173 uM by 120 min.

**FX-S1, the incumbent's 'Strecker is the only amine exit'.**

- Hofmann 2000 Table 3, within one pot: (PA + PAA) / Phe consumed 0.181 at 10 min (-0.74 dex) to 0.058 at 500 min; one-sided p under omega on each value 0.105; break-even SD 0.45 dex. Flag: the Phe column coincides with Hofmann 2000b's glucose / Phe column (dossier caveat).
- Piloty 1979 Ala loss vs Zhou's DMP under yield 1: gap 2.51 dex (pH 6.97; an upper-bound MGO exposure); p 0.046 (no buffer term) / 0.196 (phosphate at M5's kappa edge), cross-laboratory SD 1.49 dex (tau_ABS and omega, the largest listed); break-even between-laboratory SD 0.99 / 0.35 dex.
- Holm: {"piloty_vs_zhou": 0.092, "hofmann_within_pot": 0.105}; falsified: **False**.

**FX-S2, 'the amine-dependent MGO loss is a 1:1 amine-consuming adduct' (Kanzler 2026, one study):** MGO lost 95.7 mM by 120 min with Ala, at most 35.1 mM by the MGO-alone path, so >= 60.6 mM amine-dependent; a 1:1 adduct leaves <= 39.4 % Ala, observed 100 % (SD used 3); t(2 df) 20.2, p 0.0024; falsified: **True**. The Piloty-sized sink on the same course consumes 1.8 % Ala (consistent with ~0 %).

**D2, Adams 2008 pH 4 vs 7 under R1:** no contradiction: R1 needs the Strecker flux x28 at pH 4 (acetate) over pH 7 (phosphate), the buffer confounded with pH; a prediction, not a falsification (AK flux ratio 28.1).

## Job (c): the designed decider under the adopted judge (synthetic)

Sets per truth 50, scoring draws 20000, null replicates 1000 per rival, adequacy replicates 100, per-reading CV 0.08 dex. P(a) / P(b): P(credit the generating rival) with the gate passing; 95 % Clopper-Pearson in brackets; false credit: the worse direction. Gate: the M1.1 misspecified false credit k / n with the Wilson lower end (must not exceed 0.05), for the correctly specified generator (x0: a run-level reproducibility term cancels in a within-run unit) and an unmodelled per-point term of 0.5 / 1 x the declared per-point SD. A failing gate turns every credit of that pair into 'undecided'.

| pair | design | n | P(a) [CP 95] | P(b) [CP 95] | false credit (CP upper) | gate k/n (Wilson low): x0 / x0.5 / x1 | not computable | null unresolved a / b |
|---|---|---:|---|---|---|---|---:|---|
| C1_vs_C0 | M | 2 | 0.74 [0.60, 0.85] | 0.42 [0.28, 0.57] | 0.10 (0.22) | 6/100 (0.028) pass / 4/100 (0.016) pass / 4/100 (0.016) pass | 0.00 | 0.000 / 0.000 |
| C1_vs_C0 | M | 3 | 0.70 [0.55, 0.82] | 0.62 [0.47, 0.75] | 0.04 (0.14) | 2/100 (0.006) pass / 4/100 (0.016) pass / 8/100 (0.041) pass | 0.00 | 0.000 / 0.000 |
| C1_vs_C0 | D | 2 | 0.16 [0.07, 0.29] | 0.00 [0.00, 0.07] | 0.06 (0.17) | 3/100 (0.010) pass / 2/100 (0.006) pass / 3/100 (0.010) pass | 0.00 | 0.000 / 0.000 |
| C1_vs_C0 | M+D | 2 | 0.62 [0.47, 0.75] | 0.50 [0.36, 0.64] | 0.04 (0.14) | 3/100 (0.010) pass / 4/100 (0.016) pass / 5/100 (0.022) pass | 0.00 | 0.001 / 0.000 |
| C1_vs_C0 | M+D+K | 2 | 0.00 [0.00, 0.07] | 0.66 [0.51, 0.79] | 0.00 (0.07) | 0/100 (0.000) pass / 1/100 (0.002) pass / 0/100 (0.000) pass | 0.00 | 0.043 / 0.232 |
| C1_vs_C3 | M | 2 | 0.98 [0.89, 1.00] | 0.94 [0.83, 0.99] | 0.02 (0.11) | 1/100 (0.002) pass / 1/100 (0.002) pass / 2/100 (0.006) pass | 0.00 | 0.006 / 0.000 |
| C1_vs_C3 | M | 3 | 0.96 [0.86, 1.00] | 0.96 [0.86, 1.00] | 0.00 (0.07) | 0/100 (0.000) pass / 1/100 (0.002) pass / 0/100 (0.000) pass | 0.00 | 0.005 / 0.000 |
| C1_vs_C3 | D | 2 | 0.90 [0.78, 0.97] | 0.96 [0.86, 1.00] | 0.02 (0.11) | 1/100 (0.002) pass / 0/100 (0.000) pass / 0/100 (0.000) pass | 0.00 | 0.033 / 0.030 |
| C1_vs_C3 | M+D | 2 | 0.98 [0.89, 1.00] | 0.00 [0.00, 0.07] | 0.00 (0.07) | 0/100 (0.000) pass / 0/100 (0.000) pass / 0/100 (0.000) pass | 0.00 | 0.057 / 0.035 |
| C1_vs_C3 | M+D+K | 2 | 0.00 [0.00, 0.07] | 0.00 [0.00, 0.07] | 0.00 (0.07) | 0/100 (0.000) pass / 0/100 (0.000) pass / 0/100 (0.000) pass | 0.06 | 0.297 / 0.455 |
| C0_vs_C3 | M | 2 | 0.98 [0.89, 1.00] | 0.94 [0.83, 0.99] | 0.06 (0.17) | 3/100 (0.010) pass / 7/100 (0.034) pass / 5/100 (0.022) pass | 0.00 | 0.003 / 0.001 |
| C0_vs_C3 | M | 3 | 0.96 [0.86, 1.00] | 0.96 [0.86, 1.00] | 0.02 (0.11) | 2/100 (0.006) pass / 0/100 (0.000) pass / 1/100 (0.002) pass | 0.00 | 0.000 / 0.000 |
| C0_vs_C3 | D | 2 | 0.94 [0.83, 0.99] | 0.92 [0.81, 0.98] | 0.00 (0.07) | 0/100 (0.000) pass / 0/100 (0.000) pass / 0/100 (0.000) pass | 0.00 | 0.030 / 0.012 |
| C0_vs_C3 | M+D | 2 | 0.88 [0.76, 0.95] | 0.00 [0.00, 0.07] | 0.00 (0.07) | 0/100 (0.000) pass / 0/100 (0.000) pass / 0/100 (0.000) pass | 0.00 | 0.053 / 0.021 |
| C0_vs_C3 | M+D+K | 2 | 1.00 [0.93, 1.00] | 0.00 [0.00, 0.07] | 0.00 (0.07) | 0/100 (0.000) pass / 0/100 (0.000) pass / 0/100 (0.000) pass | 0.00 | 0.114 / 0.044 |
| S1_vs_S0 | plain+phosphate | 2 | 0.34 [0.21, 0.49] | 0.90 [0.78, 0.97] | 0.04 (0.14) | 2/100 (0.006) pass / 1/100 (0.002) pass / 2/100 (0.006) pass | 0.00 | 0.000 / 0.049 |
| S1_vs_S0 | plain+phosphate | 3 | 0.96 [0.86, 1.00] | 0.92 [0.81, 0.98] | 0.00 (0.07) | 0/100 (0.000) pass / 0/100 (0.000) pass / 0/100 (0.000) pass | 0.00 | 0.000 / 0.044 |
| S1_vs_S0 | plain only | 2 | 0.90 [0.78, 0.97] | 0.92 [0.81, 0.98] | 0.06 (0.17) | 6/100 (0.028) pass / 4/100 (0.016) pass / 3/100 (0.010) pass | 0.00 | 0.001 / 0.031 |
| S1_vs_S0 | phosphate only | 2 | 0.94 [0.83, 0.99] | 0.90 [0.78, 0.97] | 0.04 (0.14) | 3/100 (0.010) pass / 2/100 (0.006) pass / 2/100 (0.006) pass | 0.00 | 0.000 / 0.020 |
| S1_vs_S5 | plain+phosphate | 2 | 0.92 [0.81, 0.98] | 0.96 [0.86, 1.00] | 0.00 (0.07) | 0/100 (0.000) pass / 2/100 (0.006) pass / 2/100 (0.006) pass | 0.00 | 0.003 / 0.025 |
| S1_vs_S5 | plain+phosphate | 3 | 0.98 [0.89, 1.00] | 0.98 [0.89, 1.00] | 0.00 (0.07) | 0/100 (0.000) pass / 0/100 (0.000) pass / 2/100 (0.006) pass | 0.00 | 0.003 / 0.023 |
| S1_vs_S5 | plain only | 2 | 0.96 [0.86, 1.00] | 0.96 [0.86, 1.00] | 0.04 (0.14) | 3/100 (0.010) pass / 6/100 (0.028) pass / 5/100 (0.022) pass | 0.00 | 0.000 / 0.030 |
| S1_vs_S5 | phosphate only | 2 | 0.00 [0.00, 0.07] | 0.00 [0.00, 0.07] | 0.00 (0.07) | 0/100 (0.000) pass / 0/100 (0.000) pass / 0/100 (0.000) pass | 0.00 | 0.000 / 0.000 |
| S5_vs_S0 | plain+phosphate | 2 | 0.90 [0.78, 0.97] | 0.92 [0.81, 0.98] | 0.04 (0.14) | 4/100 (0.016) pass / 4/100 (0.016) pass / 3/100 (0.010) pass | 0.00 | 0.000 / 0.017 |
| S5_vs_S0 | plain+phosphate | 3 | 0.98 [0.89, 1.00] | 0.86 [0.73, 0.94] | 0.00 (0.07) | 0/100 (0.000) pass / 0/100 (0.000) pass / 2/100 (0.006) pass | 0.00 | 0.000 / 0.017 |
| S5_vs_S0 | plain only | 2 | 0.00 [0.00, 0.07] | 0.00 [0.00, 0.07] | 0.00 (0.07) | 0/100 (0.000) pass / 0/100 (0.000) pass / 0/100 (0.000) pass | 0.00 | 0.000 / 0.000 |
| S5_vs_S0 | phosphate only | 2 | 0.94 [0.83, 0.99] | 0.96 [0.86, 1.00] | 0.04 (0.14) | 3/100 (0.010) pass / 4/100 (0.016) pass / 0/100 (0.000) pass | 0.00 | 0.000 / 0.022 |

## Declared limits

- The rivals' structures are exact in the mini-model, so every power is an upper bound for real, misspecified data.
- MGO is held at its charge (the design measures it and forces the model with it); forced inputs are therefore not redrawn by the null (case_inputs: forced_inputs False). With a measured MGO course the null must redraw it.
- The spent rows enter only through R1's band (ks >= job (a)'s ks_min; kS near Hofmann's implied rate), not as conditioning rows of the density: a declared simplification of p(y_dec | y_spent).
- The heat-up is a measured design input (1 min, every rival alike); an unmeasured heat-up mimics R1's lag.
- The AK trap reading is (recovery x AK + blank) / reporting limit (0.5 uM), recovery 0.1 - 1 uncalibrated, blank 0.2 - 1 x the limit: a non-detect is a blank-level reading, not a censored one.
- R1's exit (ks) and the ALA -> aminoacetone transfer (k2 +- 1 dex) are declared, unsourced widths; R3 and C0 predict the step units structurally (constant-free).
- Scoring ensembles come from a defensive proposal (0.9 prior + 0.1 edge-weighted; log p / q passed to the judge); unresolved bootstrap replicates and observed sets are re-estimated once on independent ensembles 3x larger.
- The precision requirement is the shared floor (ESS >= 100), declared; the misspecified generator's extra term is a multiple of the per-point SD because the pyrazine lane has no within-laboratory term.
