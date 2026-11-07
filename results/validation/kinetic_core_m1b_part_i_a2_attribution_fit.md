# M1b part (i): A2's attribution fit F-A (prereg v3.8 sec. 21.2)

Status 2 (converged True), failed pot runs 0, 258 evaluations in 1 segment(s), 2214 s.

## Bands

- log10 k_cys_thermal (145 C, /min): [-2.818, -1.980] (Kang's total 1.048e-02 /min; rho 0.145)
- mu_h2s: [-1.048, -0.083] (Zheng pH 7 at 140 C 1.047e-02 vs Kang's total 8.648e-03 /min; class Ea min 29.7 kJ/mol)

## Result

- k_cys_thermal at fit -2.818; mu_h2s -0.221; lumped yield 3.72e-20; ARP pKa 7.32
- Meynier cysteine drifts at fit (4.5 .. 6.5): -0.500 / -0.642 / -0.548 / -0.373 / -0.357
- meets Meynier at sigma_max (0.135): **False**; by choice: sd_0.071 False, sd_0.100 False, sigma_max False, sigma_max_edge0.25 False, beta_s_x0.5 False, beta_s_x1 False, beta_s_x2 False, shape_unscaled_s_x1 False, shape_unscaled_s_x2 True
- B60 cost (52 rows + prior): 28.0863 at F-A vs 21.8833 at the confirmation fit (change +6.2030)
- band edges active: ['Ea_decay_thiol_sink', 'k_cys_thermal', 'ph_acid_yield_per_sink_event']

## Diagnostic (no verdict): r_cys_thermal and ch_cys_ox with their amine kept

- at F-A: -0.256 / -0.522 / -0.513 / -0.360 / -0.348
- at the confirmation fit: +0.069 / -0.365 / -0.490 / -0.359 / -0.349
