# Wave B58 pre-registration — TTCA formation as a measured equilibrium

*Written 2026-10-27 before any B58 number existed. Kind 3 under `scripts/generators/WAVES.md` ("How a
change is accepted"): a new reaction, so it is accepted on pre-registered hold-out tests. Generators:
`scripts/generators/generate_kinetic_core_b58_fit.py`, `scripts/generators/generate_kinetic_core_b58_ship_rule.py`.
Switch: `parameters_sulfur.TTCA_FORMATION` (`enabled`). Baseline: what ships (B51).*

## 1. Why

The sulfur lane can only BREAK TTCA (2-threityl-thiazolidine-4-carboxylic acid) down: `r_ttca_cys` (ring-
opening to cysteine + pentose, fitted) and `r_ttca_deg` (to 1-deoxypentosone, fitted). Nothing forms it. But
cysteine and a pentose condense to TTCA quickly and almost completely:

- **Cao et al. 2017** (`cao2017_extraction.md`; Cys + Xyl 0.40 M each, 0.2 M phosphate pH 5.5 → 4.15, 90 °C,
  0.5–5 h): [TTCA]/([Cys][Xyl]) = 58.5–64.6 L/mol at all six sample times — an apparent equilibrium constant
  K ≈ 60 L/mol (band 45–75), established within 30 min, with TTCA holding ~89 % of the cysteine.
- **Zhang et al. 2025b** (`zhang2025b_extraction.md`; Cys + Xyl 40 + 40 mM, 0.2 M phosphate pH 6.5, 120 °C):
  ~45 % of the cysteine in TTCA (18 mM) at 40–50 min, ~4 % as the Cys-Amadori compound (Fig. 2, digitised).
- **Zhai et al. 2019** (`zhai2019_extraction.md`; Xyl 165 + Cys 82.6 mM, water, pH 7.4 set with NaOH,
  unbuffered): TTCA plateaus at 54–63 mM (65–77 % of the cysteine) at 70–90 °C.

So every engine pot that charges cysteine with a pentose treats all its cysteine as free when most of it is
bound. That changes how much cysteine is free to release H2S, to catalyse the pentose's enolisation, and to
be lost to the sinks.

## 2. What changes

- **New reaction `r_ttca_form`: Cys + PENT → TTCA**, its rate DERIVED by detailed balance against the existing
  ring-opening, k_form(T) = K × k_ttca_cys(T), with K = 0.060 L/mmol (Cao 2017). No new fitted parameter.
- **Declared assumption:** K is held constant across temperature and pH (measured at 90 °C, pH ~4.5–5.5).
  Thiazolidine formation needs the free amine and the thiol, so K should depend on pH; nothing in the corpus
  measures how. The two mechanism tests below transfer K to 120 °C / pH 6.5 and to 90 °C / pH 7.4.
- **Candidates, fixed order of preference:** (A) refit — B51's 52-row objective with the step on, same free
  set, bands, sigmas, optimiser, budget and two-start protocol (start 0 = B51's optimum); (B) B51's frozen
  vector with the step on.

## 3. Ship rule (`generate_kinetic_core_b58_ship_rule.py`), each arm in a child process

- **T1 (mechanism, both must hold):** (a) Zhang 2025b's pot at 50 min: predicted TTCA within a factor 2 of
  45 % of the charged cysteine (22–90 %); (b) Zhai 2019's pot at 90 °C, 160 min: predicted TTCA within a
  factor 2 of 58 mM (29–116 mM, capped by the 82.6 mM charge).
- **T2:** Maillard-path hold-out sulfur rows: within-3x not lower than the baseline's by more than 2, and the
  median |log10 error| not worse by more than max(0.1 dex, one bootstrap standard error of the baseline
  median).
- **T3:** Hofmann ribose pH 5 + xylose pH 5 rows: median |log10 error| not worse by more than the same
  tolerance (these sit ~100x low under B51; recorded either way).
- **T4:** directional panel, independent claims agreeing: not fewer than the baseline's minus 3 (binomial
  noise on ~31 claims).
- **SHIP (A)** if (A) passes T1–T4, else **(B)** if (B) does, else **DO NOT SHIP** (the step stays off).

## 4. Expectations

Binding most of the cysteine lowers free cysteine several-fold in every cysteine + pentose pot, so H2S release
and the cysteine-catalysed enolisation slow, and the refit will raise them. T1(a) is the riskiest: K measured
at pH ~4.5 may be too large at pH 6.5, over-binding. Probability of SHIP: **about 0.5**.

### 4a. Found while checking the harness, before either fit start ran

A probe of B51's shipped fit with the step off and on (engine path, the T1 pots plus Cao's own pot):

| pot | measured | step off | step on |
| --- | --- | --- | --- |
| Cao 2017, 90 °C, 60 min, Cys + Xyl 400 + 400 mM | TTCA ≈ 356 mM, free Xyl ≈ 77 mM | TTCA 0, PENT **0**, Cys 155 | TTCA 1.7, PENT 0 |
| Zhang 2025b, 120 °C, 50 min | TTCA 18 mM | TTCA 0, PENT 0, Cys 0.27 | TTCA 0.0 |
| Zhai 2019, 90 °C, 160 min | TTCA ≈ 58 mM | TTCA 0, PENT 57, Cys 0 | TTCA 0.02 |

The step does what its K says (apparent K 0.054–0.059 L/mmol where pentose survives; 0.060 put in). But the
lane CONSUMES THE PENTOSE (or, in Zhai's 2:1 pot, the cysteine) long before TTCA can hold it: at 90 °C the
xylose Cao measured still present — most of it bound in TTCA — is gone in the model within the hour. The direct
cysteine-catalysed enolisation (`r_pent_dpo` / `r_pent_tdp`, PENT + Cys → DPO / TDP, fitted at 145 °C only)
probably stands in for the real route PENT + Cys ⇌ TTCA → Cys-Amadori → deoxyosone, and extrapolated down to
90 °C it is far too fast. T1 is therefore expected to FAIL unless the refit slows those steps; the probability
of SHIP in §4 drops to **about 0.15**. The finding is the next structural target either way: route the
pentose's cysteine-catalysed loss THROUGH TTCA, and test the lane's low-temperature kinetics on Cao's 90 °C
cysteine and xylose time courses (B55's single-temperature problem).

## 5. Verdict

*(Written after the run, below this line.)*

**DO NOT SHIP** (`kinetic_core_b58_ship_rule.json`, 2026-10-28). The step stays off (`TTCA_FORMATION.enabled = False`).

| test | (A) refit, cost 20.256 (B51: 21.035 start) | (B) B51 frozen, step on |
| --- | --- | --- |
| T1 Zhang 2025b TTCA, band 9–36 mM | 0.0005 mM — FAIL | 0.0005 mM — FAIL |
| T1 Zhai 2019 TTCA, band 29–82.6 mM | 0.019 mM — FAIL | 0.019 mM — FAIL |
| T2 hold-out sulfur, within 3x / median dex (baseline 1/14, 1.61) | 2/14, 0.99 — pass | 2/14, 1.07 — pass |
| T3 Hofmann pH 5, median dex (baseline 1.80) | 0.71 — pass | 0.43 — pass |
| T4 directional (baseline 17/31) | 17/30 — pass | 16/31 — pass |

As §4a foresaw, T1 fails by three orders of magnitude in both arms: the lane consumes the pentose (and in Zhai's 2:1 pot
the cysteine) in seconds, long before an equilibrium with K ≈ 60 L/mol can hold either. The refit could not slow the
trunk because no fit row measures a timescale. The informative side result is T3: even the trace of cysteine TTCA
shelters raises the Hofmann pH 5 thiols by 10–25x, which says those rows are governed by how long free cysteine
survives the first seconds.

**Caveat found the same day, before this verdict was written:** every number above (baseline included) was computed on
the 41-node uniform pH grid, which a flux attribution showed is not converged in exactly these pots (Hofmann pH 5 MFT
moves 39x when the transient is resolved). The verdict does not depend on it — T1's gap is ~1000x — but the T2/T3
magnitudes do. The structural fix (route the pentose's cysteine-catalysed loss through TTCA and the Cys-Amadori, and
constrain the trunk with measured time courses) is B59; the numerics fix is B60 (`kinetic_core_b60_prereg.md`).
