# resolve2_pyrbuf -- pre-registration (written 2026-11-01, BEFORE any test below was computed)

Scope: E5 / CF-011 (form of buffer catalysis), E4 (Strecker-step pH law vs Piloty 1979), E13 / E14 (condition
dependence of the pyrazine mixing ratio R = MPZ^2 / (PZ x DMP)), early steps (Huber F15 -> F07; the a_w maximum; Wolf's
lysine plateau, F05 / CF-072). Split: docs/reference/contrast_benchmark/split_v2_per_change.yml. TEST contrasts used:
E5 {potman1989, kallen1971b, thornalley1999} + the keyed TEST table heldout_larralde1995.yml; E13 {lee2024 keyed TEST
table, shibamoto1977, low2007, vanlancker2012, amranihemaimi1995}; E4: no TEST lab under v2 except F16#11 (monforte2021);
Piloty 1979 is a NEW laboratory (TU Berlin, Baltes) not in the frozen split, so it is unsplit, not TEST. Bands come only
from DEVELOPMENT labs or from unsplit papers, never from the TEST rows they are scored on.

Data classes for E5 (from the split): DEV labs = auburn_bell (bell1995/1997), baynes_south_carolina (zyzak1995),
ghent, leuven, nestle (davidek2002), reading, rizzi_miami (rizzi2004/2010), scut_guan_yu (guan2012), wageningen.
UNSPLIT (new since the freeze): saunders1966 (Cerebos), davies2005 (Delaware / Leeds lineage). TEST: potman1989,
kallen1971b, thornalley1999, larralde1995.

---------------------------------------------------------------------------------------------------------------------
## Q1. The form of buffer catalysis (E5, CF-011)

### Rivals (stated neutrally, with chemistry)
Form of the concentration law, per step:
- **F-lin** -- general acid / base catalysis of one elementary step: k = k_u(pH, T) + sum_b kcat_b c_b, linear in each
  buffer's total concentration (Bronsted; proton transfer to / from the enol / enediol in the catalysed TS). Written as
  k = k_u (1 + kappa_b c_b), kappa_b in 1/M.
- **F-sat** -- the E5 fix as written: a Michaelis-type factor k = k_u (1 + A c / (K_B + c)) with half-saturation inside
  the food range (K_B <= 0.1 M), e.g. buffer binding to the substrate or catalyst site saturation.
- **F-net** -- linear per elementary step (F-lin), with APPARENT saturation only where a preceding uncatalysed step
  becomes rate-limiting (change of rate-determining step, Jencks / Kallen). The ceiling is then set by the uncatalysed
  step and is predictable from data that do not involve the buffer.
- **F-null** -- no buffer-species catalysis; apparent effects are pH control, trace metals carried by the salt, or
  ionic strength.
Anion specificity:
- **S-P** -- phosphate-specific (bifunctional H2PO4- / HPO4 2-, or an XO2- motif): carboxylates and bicarbonate ~0.
- **S-gen** -- general acid / base: every buffer acts, scaled by its active-form concentration and pKa (Bronsted).
- **S-step** -- anion selectivity differs by STEP: carboxylates act on some steps (free-sugar enolisation) and not on
  others (Amadori decay).
pH: **pH-ratio** -- the catalysed increment kcat c is pH-flat or broad and the RATIO k/k_u fades only where k_u rises
(base-catalysed k_u at pH >= 8); **pH-species** -- the increment tracks one phosphate species fraction.
Browning vs Amadori decay: **B-same** (the browning factor IS the Amadori-decay factor) vs **B-chain** (browning
integrates >= 2 catalysed steps, at least one not the Amadori decay).

### Bands (non-scored data only)
kappa_P (per M TOTAL phosphate) from DEV data, any amine, 25 - 130 C, pH 5 - 7.4: Guan 2012 DFP decay 110 C x3.8 at
0.1 M (0-10 min) / x4.4-4.8 (60-min mean); 130 C x1.5 / x2.6; Davidek 2002 DFG decay 90 C pH-stat 0.1 M x3.0 / 3.5 /
2.06 / 1.27 at pH 5 / 6 / 7 / 8; Zyzak 1995 (37 C, chelator, N2) kappa 15; Bell 1997 (25 C) 7 - 17. => DEV band
kappa_P in [5, 40] /M at pH 5 - 7.4 (central ~20); pH 8 [2, 10] (Davidek only). Carboxylate on Amadori decay: Guan
acetate / citrate = water within +-15 % at 0.1 M -> kappa_carboxylate(ARP decay) in [0, 1.5] /M.

### Tests and decision rules (fixed now)
- **T1 Potman (TEST), form.** Net catalysed glucose-consumption increment (5-h secant, Fig. 4; and the printed log-log
  slope 1.34) between 0.1 and 0.5 M phosphate. F-sat (K_B <= 0.1 M) predicts increment ratio I(0.5)/I(0.1) <= 1.67;
  F-lin predicts 5. Rule: observed ratio >= 3 (after propagating +-1.5 % read error) rejects F-sat. Power check first:
  synthetic Potman data from F-sat (K_B = 0.03 and 0.1 M) and F-lin at the read noise; the rule must classify >= 90 %
  correctly, else T1 is not scored.
- **T2 Potman and Larralde (TEST), size.** kappa from T1 (Potman, 100 C, pH 5.6, sugar consumption with Gly) and
  kP / k0 from Larralde's printed rate law (100 C, pD 7.4, amine-free ribose). PASS if inside the DEV band [5, 40].
- **T3 Larralde (TEST), specificity.** Bicarbonate / phosphate at 0.05 M = 0.52. S-P predicts k0 / k(0.05 M phosphate)
  = 0.38 (k0 band x[0.79, 1.26] for a 0.1 dex extrapolation error -> 0.30 - 0.48). Observed above 0.48 -> bicarbonate
  catalyses (S-P rejected for this step); report bicarbonate's per-mole activity relative to phosphate.
- **T4 Kallen (TEST), saturation.** F-net predicts the phosphate ceiling (pH 6.1 - 6.6) = (attack-limited rate per
  amine free base, read from Kallen Fig. 6's buffer-free pH 4.5 - 5.2 plateau) / (buffer-free rate at pH 6.1 - 6.6).
  PASS if within 0.3 dex of the observed x3.5 - 5. F-sat has no independent prediction of the ceiling (it fits it).
- **T5 Thornalley (TEST), shape only.** F-sat predicts a concave dose response 10 -> 100 mM; observed convex or linear
  rejects F-sat. kappa from Thornalley is NOT scored (metal confound, never crossed with the chelator; 37 C).
- **T6 pH (DEV Davidek + TEST Potman), diagnostic.** Additive form predicts ratio(9) for Potman = 1 + I(7) / k_u(9)
  from its own pH-7 increment; PASS within 0.15 dex of the printed x3.6 means the "fade at pH 8 - 9" is the k_u rise,
  not a loss of catalysis.
- **Negative control.** Guan's acetate / citrate (DEV) on DFP decay and Bell's citrate at 25 C must come out kappa ~0
  under the same arithmetic; Rizzi's bis-tris is the inert reference. If the arithmetic gives the carboxylates
  kappa > 5 on Amadori decay, the method is broken.
- **B-same vs B-chain (all labs).** B-same predicts that carboxylates catalyse browning only where they catalyse
  Amadori decay; B-chain allows browning catalysis without Amadori-decay catalysis. Also B-same predicts browning order
  ~1 in buffer; orders > 1.2 (Bell 1.6, Rizzi 2.5 on A420, Potman 1.34) favour B-chain.

### Expectations (before computing)
| claim | P |
|---|---|
| T1 rejects F-sat (power adequate) | 0.80 |
| T2 Potman kappa in [5, 40] | 0.70 |
| T2 Larralde kP / k0 in [5, 40] | 0.65 |
| T3 bicarbonate / phosphate above the S-P interval | 0.60 |
| T4 Kallen ceiling within 0.3 dex of F-net | 0.55 |
| T5 Thornalley not concave | 0.75 |
| T6 Potman fade reproduced by the additive form | 0.60 |
| Verdict on form: F-lin / F-net (linear per step, network saturation) | 0.70 |
| Verdict on specificity: S-step | 0.60; S-gen 0.30; S-P 0.10 |
| Browning is NOT just the Amadori-decay step (B-chain) | 0.75 |

---------------------------------------------------------------------------------------------------------------------
## Q2. Strecker-step pH law (E4) vs Piloty 1979 (owner decision stands: no change; declared band [0, 0.58])

Rivals: **E4-step** (shipped: the below-7 slope 0.58 dec / pH sits on the Strecker step k_mgo_ak / k_go_ak) vs
**E4-supply** (the step is flat, s ~ 0-0.15; the cascade's slope belongs to the dicarbonyl supply / condensation).
Piloty feeds MGO (no sugar), so the supply is removed: E4-supply predicts flat consumption, E4-step predicts the
Strecker part of consumption x0.038 from pH 6.97 to 4.52.
- Auxiliary that decides relevance: the Strecker share f of measured amino-acid consumption. The test discriminates
  only if f is not small. **Test Q2-a:** bound f from Piloty's own pH-flatness under E4-step (k_tot(4.52)/k_tot(6.97)
  = 1 - f (1 - 0.038), +-3 % loss error). **Test Q2-b (engine, negative control included):** the engine's glycine loss
  in Piloty's pot (MGO 312 mM, Gly 16 mM, 100 C, 3 h, clamped pH 4.52 / 6.97 / 9.5) as a share of the measured loss;
  control: no MGO -> no glycine loss. If the engine's step is <= 1 % of the measured loss, Piloty constrains a channel
  the engine does not carry, and is NOT discriminating for E4 (the step law can hide under any non-Strecker adduct
  channel).
- Expectation: Piloty discriminates step vs supply only if f >= 0.3: P 0.25. Expected verdict: evidence in the E4-supply
  direction only under an auxiliary (pH-independent Strecker share), weight low; record for a future revision.

---------------------------------------------------------------------------------------------------------------------
## Q3. Condition-dependent R (E13 / E14; P1x / P1y / P2b)

Rivals:
- **P0** one shared aminoketone pool and one condensation constant with statistical factors: R = 4 everywhere.
- **P1x / P1y** pair-specific constants (k_AA, k_AB, k_MM): R = 4 k_AB^2 / (k_AA k_MM) = 4 r^2. Condition dependence
  needs the pairs to differ in pH / T / N dependence (different Ea or pKa sensitivity per pair).
- **P2b** pairing at R_pair = 4 plus an extra MPZ route with its own condition dependence: MPZ = MPZ_pair + E. Three
  sub-forms by what E needs: E-A (aminoacetaldehyde side, E ~ a), E-M (aminoacetone side, E ~ b), E-0 (pool-independent).
- Operator rivals: co-elution of 2,5- / 2,6-DMP (Ghent, Koehler), vapour partition (Lee; R is invariant to a
  log-linear volatility series PZ -> MPZ -> DMP), the Lee MPZ axis x10 (text vs figure).

The pattern that separates them (derivation, before any number): with pools a (aminoacetaldehyde) and b (aminoacetone)
and fixed constants, PZ ~ a^2, DMP ~ b^2, MPZ_pair ~ 2 r a b, so under P1 R = 4 r^2 does NOT depend on a / b (pool
invariance). Under P2b, R depends on the pools: E-A gives R = (2 + e / b)^2, E-M gives R = (2 + e / a)^2, E-0 gives
R = (2 + E / ab)^2. **Test Q3-a (supply perturbation at fixed chemistry):** within a pot's time course, or an arm that
changes the supply without changing pH / T / N source, calibrate P2b's extra share from the first point's R and predict
R at the other points from the measured PZ and DMP (a ~ sqrt PZ, b ~ sqrt DMP); P1 predicts R unchanged. Data: Koehler
1969 Asn +- acetaldehyde (unsplit); Lee 2024 200 C time course 22 / 35 / 54 / 88 min and 170 C 54 / 88 min (TEST; the
x10 MPZ issue cancels in R ratios); Shibamoto 1976 100 C time course 30 min - 48 h (unsplit). Rule: a rival whose
predicted R ratio lies outside the observed ratio's band (read error / SD propagated, 2 sigma) at a point fails it.
Condition series (Bemis-Young pH, Wong N source, Shibamoto T) change the chemistry, so they cannot separate P1 from P2b
(both have free condition laws); they are scored only for P0 (R = 4 within 0.3 dex).
- **Test Q3-b (Lee level):** P0 / pairing floor R <= 4 h (h <= 1) at face value and with MPZ / 10.
Expectations: P0 rejected by >= 3 labs (P 0.95); Q3-a separates P1 from at least one P2b sub-form (P 0.5); Lee's level
verdict flips with the /10 sensitivity (P 0.8).

---------------------------------------------------------------------------------------------------------------------
## Q4. Early steps

- **Huber -> F07.** Rival H-redox (the 1-amino-1,4-dideoxy-2,3-hexodiulose is a reduction product, F15) vs H-elim
  (C-4 beta-elimination from the 2,3-enediol, amine retained: F07). Test: oxidation level of the structure (same as the
  1-deoxyosone -> elimination). Expect H-elim P 0.85; split the McGill cell (the RA[4,5] share -> F07; pyruvaldehyde ->
  1,2-propanediol stays F15).
- **a_w maximum (Cui 2026 vs Pereyra / Warmbier / Eichner).** Rivals: W-max (a rate maximum at a_w 0.5 - 0.75 in
  every matrix) vs W-matrix (the location and existence of a maximum depend on matrix: humectant, glass, water content;
  Pereyra's powder is flat 0.31 - 0.71) vs W-endpoint (Cui's 20-min end points of consumed intermediates, a_w not held,
  cannot show a rate maximum). Test: does any pair measure the same channel (rate of lysine / amine loss) with a_w held?
  Expect W-endpoint for Cui (P 0.85): no engine change.
- **Wolf plateau (F05, CF-072).** Rivals: L-noreturn (Nguyen: lysine consumed only by glycation; no plateau with sugar
  in excess); L-return (the engine's r_flp_decay writes lysine back: a pseudo-steady-state free-lysine fraction
  k_d / (k_d + k_g [Glc])); L-pool (an inaccessible lysine fraction); L-sugar (glucose depleted); L-assay (FDNB counts part
  of the bound Amadori lysine). Test: second-pass prediction of each (Wolf Table 2: no further loss), with the frozen
  B20 constants at 125 C. Expect L-noreturn fails (P 0.7).
