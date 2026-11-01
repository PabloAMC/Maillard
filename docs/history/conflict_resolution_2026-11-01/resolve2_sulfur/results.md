# resolve2_sulfur results (2026-11-01)

Prereg: `prereg.md` (secs. 1-4 before any check; Amendment 1 compute, before results; Amendment 2 POST HOC).
Numbers: `phi1_retest.json`, `phi1_face.json` (post hoc). Code: `phi1_retest.py` (imports the cluster-5 module
unchanged, swaps bands only), `phi1_face.py`. Text proposals: `corrections.md`; register: `register_proposals.yml`.

## 1. CF-076 / E10 (DEV-lab feasibility only; E10 has one TEST lab, Mizutani 1960, UNTESTABLE on held-out)

### 1a. The release band at 100 C, from data other than T1-T8

k_r(100 C, pH 4.5) = k25_6 x 10^(-1.1 s_lo) x exp(Ea/R (1/298.15 - 1/373.15)).
- k25_6 (uncatalysed release of a pentose-Cys ring, 25 C, pH 6.0): Karpinski 1990 AThz 1.8e-4 /min (RSD 30 %); AThz /
  Thz = 1.0 (phi ~1, Warsaw); Karpinski's TC is x3.3 Kallen's corrected TC (5.4e-5), read as phosphate / ionic-strength
  inflation of his pH 6 runs -> declared band [3e-5, 2.5e-4].
- s_lo (acid limb, 4.5 -> 5.6): Kallen TC 0.54 +-0.3; Karpinski PrThz 0.49; Bustillo 2018 (aryl / alkyl rings, 20 C)
  puts the maximum at 4.4-5.2, i.e. s_lo ~ -0.2 to -0.6 (makes k_r45 LARGER).
- pKa_T: Radomski 1989 charge balance 5.2-6.1 (no measurement).

| Ea (kJ/mol) | x(25->100 C) | k_r45(100), face value (k25 1.8e-4, s_lo 0.49) | full band (k25 3e-5..2.5e-4, s_lo 0.85..0) |
|---|---|---|---|
| 20 | 5.1 | 2.6e-4 | 1.8e-5 - 1.3e-3 |
| 30 | 11.4 | 5.9e-4 | 4.0e-5 - 2.8e-3 |
| 40 | 25.6 | 1.3e-3 | 8.9e-5 - 6.4e-3 |
| 50 | 57.6 | 3.0e-3 | 2.0e-4 - 1.4e-2 |
| 64 | 179 | 9.4e-3 | 6.2e-4 - 4.5e-2 |
| 80 | 656 | 3.4e-2 | 2.3e-3 - 1.6e-1 |

The cluster-5 declared band 2-8e-4 /min needs Ea 17-34 kJ/mol at face value; Ea <= ~50 at the band's slowest corner
(k25 3e-5, s_lo 0.85); never Ea 64 unless k25 < 3e-5.

### 1b. Re-test (S1 = T1, T2, T4, T5, T6, T8a, T8b; DE maxiter 20, popsize 12, seed 3)

| run | best max abs z | where it sits | reading |
|---|---|---|---|
| NC1 / NC2 (20 draws) | 58.2 / 33.7 | -- | controls fail as designed |
| positive control: E4c UNcoupled, pKa_T 5.2-6.1 | 1.29 | pKa_T 5.21 (new floor), Ea 43 | wrapper and maxiter 20 fine (cluster 5: 1.45) |
| **E4c coupled (phi ~1), Ea 40-100** | **1.76** | **k25_6 3.05e-5 (band floor), Ea 40.6 (floor)**, s_lo 0.38, pKa_T 5.55 | feasible only at the double corner |
| **E2c coupled, Ea 40-100** | **1.40** | **k25_6 3.2e-5 (floor), Ea 42.9**, s_lo 0.80 (near ceiling 0.85), pKa_T 5.38 | feasible only at the corner |
| E4c coupled, Ea 64 fixed | 11.95 | T4 11.3, T8b 12.0 | REJECTED |
| E4c coupled, Ea 30 fixed | 1.10 | k25_6 3.1e-5 (floor), s_lo 0.55, pKa_T 5.29 | feasible |
| E4c coupled, Ea 50 fixed | 3.43 | k25_6 3.0e-5 (floor); T8b 3.4, T5 3.1 | REJECTED |
| POST HOC E2c, k25_6 = 1.8e-4 (printed), Ea 10-100 | 1.39 | Ea 26.5, s_lo 0.81, pKa_T 5.44 | feasible only far below class Ea |
| POST HOC E2c, k25_6 = 8.7e-5 (mid-band), Ea 10-100 | 1.29 | Ea 31.9, s_lo 0.85 (ceiling), pKa_T 5.47 | feasible only at Ea ~32 |

### 1c. Verdict (pre-registered rule)

- Coupled family feasible in the class band, but only at Ea 40-43 (its floor) WITH k25 at the band floor (Karpinski's
  printed rate deflated x6); rejected at Ea 50 (3.4) and 64 (12.0). Required Ea by anchor: 26 (printed 1.8e-4), 32
  (8.7e-5), ~41 (3e-5). Every feasible point puts pKa_T at 5.2-5.5 (inside Radomski's band): the pKa no longer binds.
- **CF-076: EXPLAINED (provisional) CONDITIONAL on a release Ea <= ~45 kJ/mol; P 0.45 -> 0.35.** Not OPEN by the rule,
  but the admission now rests on an unmeasured, low-for-class Ea and on the steep acid limb (s_lo 0.8-0.85 in 3 of 5
  feasible points). Bustillo's acid-side maximum (s_lo < 0) would make it harder still (not modelled; lenient case).
- Expectations: positive control (0.9) yes; NC (0.95) yes; feasible in class band (0.45) yes, at the floor; Ea 64
  infeasible (0.7) yes; face value infeasible (0.75, post hoc) WRONG once Ea may drop to 26.
- Caveats: maxiter 20 (positive control reached 1.29, below cluster 5's 1.45, so search power is adequate); DEV labs
  only; T3 still fails everyone.
- Decider: one Arrhenius series (25-90 C, pH 4.5 / 6) of cysteine release from TTCA or AThz.

## 2-4. See corrections.md (B: van den Ouweland; A: Wei; C: CF-031 / E7) and register_proposals.yml.
