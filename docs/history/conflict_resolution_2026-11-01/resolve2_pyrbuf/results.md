# resolve2_pyrbuf -- results (2026-11-01). Prereg: prereg.md (written before any number below).

Scripts (stdlib, run in the maillard env via ./scripts/docker_maillard.sh): `tests.py` -> `results_tests.json`;
`shib_abs.py` -> `shib_abs.json`; `piloty_engine_probe.py` -> `piloty_engine_probe.json` (engine, in memory, clamped pH;
nothing written to results/). Nothing tracked was written.

## Q1. The form of buffer catalysis (E5, CF-011)

| test (prereg) | data class | result | verdict |
|---|---|---|---|
| T1 Potman form: I(0.5) / I(0.1) | TEST | 5-h secant **5.7** (read-error range 4.6 - 7.3); 3-h 6.2 (5.0 - 8.1); Fig. 5 slope refit 1.37 (printed 1.34). F-sat (K_B <= 0.1 M) allows <= 1.67 | **F-sat rejected.** Power: the rule calls F-lin "not saturating" 100 %, F-sat (K 0.03 / 0.1 M) 0 % of 4000 synthetic designs |
| T2 size vs DEV band [5, 40] /M | TEST | Potman kappa 22 - 29 /M at every c (5 h); Larralde kP / k0 = **32.7 /M** | both PASS (weak test: the band spans 0.9 dex) |
| T3 Larralde bicarbonate / phosphate | TEST | S-P predicts 0.38 (0.33 - 0.43); observed **0.52** -> bicarbonate kcat 0.021 /M/min, **0.23x phosphate per mole** | S-P rejected for free-ribose loss (temperature of the bicarbonate run inferred, P 0.8) |
| T4 Kallen ceiling under F-net | TEST | buffer-free attack plateau 10^3.4 predicts ceilings x3.6 / 3.9 / 6.8 at pH 6.12 / 6.31 / 6.62; observed x5 / 3.5 / 4.3: **+0.14 / -0.04 / -0.20 dex** | **F-net PASS** (all within 0.3 dex). Kallen's saturation is a change of rate-determining step, not buffer binding |
| T5 Thornalley shape | TEST | GO, 3-DG, MGO all CONVEX 10 -> 100 mM (50 mM point 15 - 28 % below the chord) | F-sat rejected; kappa not scored (metal confound) |
| T6 Potman pH 9 fade, additive form | TEST + DEV | predicted ratio(9) = 1 + I(7) / k_u(9) = 4.0 vs printed 3.7 (**-0.03 dex**) | PASS: the "fade" is k_u rising, not lost catalysis |
| Negative control (carboxylate on Amadori decay, Guan 110 C) | DEV | acetate 1.6, citrate 3.4 /M (< 5) | sound |

**DEV band provenance** (kappa_P per M total phosphate): Guan DFP 110 C 28 / 38, 130 C 5 / 16; Davidek DFG 90 C pH 5 / 6 / 7 /
8: 20 / 25 / 11 / 2.7; Zyzak 37 C 15; Bell 25 C ~10. **kappa is nearly temperature-flat from 25 to 110 C** (10 -> 38):
catalysed and uncatalysed steps carry similar barriers.

**Answers.**
1. **Form: linear per elementary step (F-lin), with apparent saturation only where the network changes its rate-limiting
   step (F-net). P 0.8.** Linear: Potman (TEST, power 1.0), Davies (acetate order 1.15, unsplit), Zyzak (DEV, linear to
   0.2 M), Larralde (authors' linear law), Thornalley (convex), Saunders (none seen to 0.5 M). The one saturation on disk
   (Kallen) is predicted from buffer-free data within 0.2 dex.
2. **Specificity: step-specific (S-step) P 0.6; general (S-gen) 0.3; phosphate-only (S-P) 0.1.** Carboxylates and
   bicarbonate act on free-sugar degradation and browning (Saunders pH 3 citrate = phosphate per mole, acetate 1/3;
   Rizzi citrate 0.45, acetate 0.12; Davies acetate; Larralde bicarbonate 0.23). They do NOT act on Amadori decay (Guan:
   acetate / citrate = water at 110 - 130 C, <= 0.1x phosphate) or on glycine loss at 25 C (Bell). The Amadori-decay null
   is single-lab (Guan, proline ARP): provisional.
3. **pH.** The catalysed increment has a broad maximum at pH 6 - 7 (Davidek pH-stat 0.10 / 0.40 / 0.55 / 0.29 per h at
   pH 5 / 6 / 7 / 8; Potman 21 / 94 / 126 / 114 % per h at 3 / 5.6 / 7 / 9). The RATIO fades at pH 8 - 9 because k_u
   rises. A species law is not identified: "k_u x [H2PO4-]" matches Davidek within 0.2 dex but not Potman's pH 9; no
   single fraction (H2PO4-, HPO4 2-, their product) fits both. Carry pH as a band. Saunders' x6 - 13 at pH 7.5 - 8 is
   browning, not the Amadori step: no contradiction with Davidek's x1.27.
4. **Browning is not the Amadori-decay step (B-chain), P 0.75:** carboxylates catalyse browning without catalysing
   Amadori decay, and browning orders in buffer exceed 1 (Bell 1.6, Rizzi 2.5 on A420, Potman 1.34, Davies 1.15 - 1.3).
   Browning factors (x2 - 13) must not be transferred to single steps.

**Minimal E5 structure (proposal; A3; red team before code):**
`k_s = k_s,u(pH, T) * (1 + sum_b kappa_{s,b} c_b)` per step s, c_b = total buffer molarity, no Michaelis term; saturation
only where the network has the uncatalysed preceding step (the thiazolidine attack). Constants measured in buffer are
divided by their own (1 + kappa c_cal) before use elsewhere.
- Amadori formation (k_schiff) and Amadori decay (k_ama_tdg / _odg / _mgo, k_flp_decay): phosphate kappa [5, 40] /M at
  pH 5 - 7.4, [2, 10] at pH 8; carboxylates / bicarbonate [0, 3] /M.
- Free-sugar enolisation / amine-free degradation: phosphate [5, 40] /M; other buffers rho x kappa_P with rho in
  [0.1, 1].
- Oxidative steps (glucosone / GO, F08): no buffer factor; a metal term instead (Thornalley's chelator effect).
- T: kappa T-flat over 25 - 130 C (within the band). Interim (no code): a declared transfer band x[1/(1 + 40 c_cal),
  1/(1 + 5 c_cal)] on phosphate-calibrated constants used in water or citrate pots.
- TEST status: Potman, Larralde, Kallen and Thornalley pass. Only the form test has power; the size band is wide.

**CF-011.** At 0.2 M the band gives x2 - 9, which covers the 3 - 6x residual. That shows phosphate catalysis is
SUFFICIENT; it does not show it occurs. Kallen (no buffer catalysis of the attack below pH 5) places any catalysis in
Cao's pot on the Amadori rearrangement, not the ring step. The S-step reading adds a cheap, sharp decider: a citrate
arm at equal pH in buffer arm C. Citrate = water favours step-specific phosphate catalysis; citrate ~0.3 - 1x
phosphate favours general catalysis. P(phosphate catalysis of formation) stays 0.5; status open.

## Q2. Strecker-step pH law (E4) vs Piloty 1979

- Q2-a: under the shipped step law (x0.038 from pH 6.97 to 4.52), Piloty's pH-flat consumption allows a Strecker share
  f of the measured loss of at most **0.04 (Ala), 0.00 (Leu), 0.20 - 0.29 (Phe, Gly, Ile)** (2-sigma, +-3 % points).
- Q2-b (engine): the engine's glycine loss in Piloty's pot is 1.9e-5 / 5.0e-4 / 1.6e-3 (shipped) and 4.0e-4 / 4.0e-4 /
  1.2e-3 (s = 0) at pH 4.52 / 6.97 / 9.5. That is **1e-5 to 1e-3 of the measured 85 / 86 / 82 %**. Negative control (no
  MGO): 4.6e-4, the same floor. MGO falls 312 -> 5.2 mM through r_mgo_mel. Both rivals predict essentially zero
  glycine loss.
- **Verdict: Piloty does NOT discriminate step vs supply (P 0.85).** Its signal is non-Strecker MGO-amino acid
  consumption, a channel the engine does not carry. The step law hides under any f <= 0.04, and Hofmann's fed-MGO
  aldehyde + acid yields (~1 - 3 mol %) make a small f likely. **Evidence for a future revision:** with a pH-independent
  Strecker share (an assumption, untested), Piloty is a second laboratory, after Hofmann, for "fed-dicarbonyl amine
  chemistry is not suppressed in acid at 100 C, and is acid-favoured at 150 C". Likelihood ratio for E4-supply ~1.5;
  the owner's decision stands. The designed E4 pot (sec. 3.1) must measure the Strecker aldehyde, not amino-acid loss.
- **New engine gap** (registered below): consumption by non-Strecker dicarbonyl-amine adducts runs ~1e3x the engine's
  Strecker consumption in Piloty's pot (implied k >= 3.5e-5 L/(mmol min) at 100 C, against k_mgo_ak 3e-8). This matches
  the Sedehi 2013 analogue (~1e3x for glyoxal + glycine). The engine's r_mgo_mel consumes no amine.

## Q3. Condition dependence of R (E13 / E14)

- **P0 (one constant, R = 4)** fails every condition series. Bemis-Young: 0.020 - 3.3 over pH 6 - 9.6. Wong: 0.28 - 556
  by N source. Shibamoto: 0.31 - 3209 by T. Koehler: 0.07 - 50. Lee (TEST): 15 - 26 at face value, 0.15 - 0.26 with
  MPZ / 10.
- **Lee 2024 (TEST) breaks the pairing bounds at both readings, in opposite directions.** Face value R 15 - 26 exceeds the
  P1 bound R <= 4. With MPZ / 10, R 0.15 - 0.26 falls below CF-087's floor R >= 0.4. The x10 conflict decides the
  direction, not the failure. In R RATIOS (pool-invariance test) the x10 cancels.
- **The separating pattern is pool invariance.** With fixed constants, P1x / P1y give R = 4 r^2, independent of the
  aminoketone pools. A P2b extra route of LOWER order in the pools than the pairing (E-A ~ a, E-M ~ b, E-0
  pool-independent) makes R fall as the pools grow. Scores (calibrated on each series' first point):
  - **Shibamoto 1976, 100 C, 8 M NH3, 30 min - 48 h** (absolute = weight % x yield; yield moves x7.6): observed
    dlog R +0.00 to +0.13. E-A, E-M and E-0 predict -0.36 to -1.20 dex: |z| 3 - 6 at an assumed sigma 0.21 dex. P1 is
    exact. Exception: the 15-min point sits 0.58 dex above (low yield 0.013 %; heat-up).
  - **Koehler 1969, Asn +- acetaldehyde** (PZ / DMP x186): observed R x1.96 (band 0.88 - 4.3). E-A predicts x6.9 (fail),
    E-M x0.21 (fail), E-0 x0.71 (marginal fail), P1 x1.0 (pass). Caveat: acetaldehyde also changes the chemistry.
  - **Lee 2024, 200 C 22 -> 35 / 54 / 88 min, 170 C 54 -> 88** (TEST): sigma(log R) 0.19 - 0.65 dex; every rival
    |z| < 1.5. **Not discriminating.** The 88-min drop of 1 dex (PZ up x3.4) is 1.3 sigma.
- **Verdict:** P0 rejected (5 labs). P2b as a lower-order extra MPZ route is rejected (P 0.75; Shibamoto's time course,
  one lab, with Koehler as support). What survives is condition-dependent pair constants (P1x / P1y), or a second-order
  extra route, which is kinetically identical to P1x. Only labels (Zhou 2024 raw ions; EX-052 stays open) can separate
  those two. The constraints the survivors must meet:
  - an effective 2E_AB - E_AA - E_MM of ~160 kJ/mol from 25 to 60 C (Shibamoto, NH3), but only ~24 kJ/mol from 170 to
    200 C (Lee, Leu; ~0.2 dex);
  - ~0.6 dex per pH unit (Bemis-Young, Gly);
  - N-source order not monotone in basicity (Wong NH4OH 14.6 < NH4 formate 556).
  One shared constant with an Arrhenius law cannot do this. **For the engine: no structure change.** R gets a declared
  condition band, not a constant.
- **CF-087:** low R in glycine pots is reproduced by a second laboratory (Bemis-Young 0.02 - 0.09 at pH 6 - 7; Zhou SI
  0.34). The "low-level data" explanation drops (0.4 -> 0.2); pH / segregated formation rises.

## Q4. Early steps

- **Huber -> F07: RESOLVED (P 0.85).** The 1-amino-1,4-dideoxy-2,3-hexodiulose has the oxidation level of the 1-deoxyosone.
  It is a C-4 beta-elimination from the 2,3-enediol with the amine kept, not a reduction. Split the F15 McGill cell: its
  RA[4,5] ~5 % share moves to F07; the pyruvaldehyde -> 1,2-propanediol side experiment (redox) stays in F15. Add a
  Huber cell to F07 (operating, maltose ARP; glucose ARP tentative).
- **a_w maximum: EXPLAINED, no conflict on one channel (P 0.85).**
  - Cui 2026 gives 20-min end points. CML / CEL are consumed intermediates. The a_w is set before heating, in open
    dishes, with no pH. An end-point maximum (CML / CEL at 0.75, 150 C) is not a rate maximum.
  - Pereyra (rate, powder, 50 - 60 C) is flat over 0.31 - 0.71. Warmbier / Eichner (glycerol, syrup, 25 - 50 C) peak at
    0.41 - 0.55: the location depends on the matrix (humectant, glass, water content; Eichner).
  - Cui's direction (0.75 > 0.97, x1.4 - 2.2) agrees with the engine table's fall above 0.85. Its 0.53 < 0.75 contrast
    is outside the table's shape (the table is flat 0.43 - 0.69), at low weight.
  - Engine note (F05 sec. (d)): protein glycation has no a_w term at all. That gap is real and goes to designed data:
    protein + sugar, 120 - 150 C, a_w held, rate series.
- **Wolf plateau vs the engine: OPEN, and NOT a CF-072 decider.** With the frozen B20 constants at 125 C, matched to pass 1
  (9.6 -> 5.0 mg/g in 2 min), pass 2 loses **36 %** with the shipped lysine-returning r_flp_decay. It loses **48 %**
  with Nguyen's no-return structure. Wolf finds 0 - 9 % (n.s.). Both structures fail. The shipped return structure's
  pseudo-steady free fraction is 0.21, not the 0.52 observed. Glucose stays in 3x excess (L-sugar, P 0.1). Left: an
  inaccessible lysine fraction of ~50 % at a_w 0.65 (L-pool, P 0.5) or FDNB counting part of the bound Amadori lysine
  (L-assay, P 0.35). Decider: Wolf 1977a (constants, glucose levels) and Wolf 1978. Under L-pool the plateau depth does
  not move with glucose; under L-assay it moves with FLP.
