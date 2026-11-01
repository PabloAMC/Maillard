# Red team: clusters 2, 3, 4 (2026-11-01)

Read-only. I re-ran these through `./scripts/docker_maillard.sh run` (PYTHONPATH=worktree): cluster3 `ph_controls_inventory.py`,
`buffer_capacity.py`, `ph_arp_pka_sensitivity.py` and `cf065_speciation.py`. Their outputs match the logs on disk.
New probe: `scratchpad/redteam_c234/meynier_probe.py`. Grosch PDF text: `scratchpad/redteam_c234/grosch1992.txt`.

## Claim 1 (cluster 3): the pH model fails two controls

- **Shipped configuration: yes.** The controls call `integrate_sulfur` the way `engine.py:2597` does, with
  `core_parameters(SULFUR)` and `core_ph_drift()` (acid_yield 7.5e-5, ARP pKa 7.07). Settings: B51 all on,
  PH_COUPLING fixed_point_nodes, TRUNK_SOLUTES_IN_BALANCE off.
- **Compositions: correct.**
  - Meynier: 35 / 50 mM, the dossier's halving inference (Flag 4). Without the halving the acid load doubles, which is
    worse for the model.
  - Hou: GX-ARP 5 mM, water, 7.0, 120 C.
  - Caveat for Hou: the ARP was eluted with 0.2 M NH4OH, its purity is not stated, and the initial SID is unknown.
- **Probe (Meynier, shipped).** Cooled end pH and drift by label:

  | label | shipped | CBX set to 0 | amine carried as NH4+ (pKa 9.25) |
  |---|---|---|---|
  | 4.5 | 3.71 (-0.79) | 0.00 | +0.30 |
  | 5.5 | 4.85 (-0.65) | 0.00 | +0.05 |
  | 6.5 | 6.11 (-0.39) | 0.00 | +0.03 |

  - CBX at the end is 21.9 / 37.6 / 42.2 mM. ACID x yield is about 2 uM, so CBX is the model's ONLY acid source.
  - The model's 0.2 M pyrophosphate (GKL + Davies) tolerates only 3.9 / 13.8 / 20.8 mM of formic / acetic acid before
    the drift passes 0.2. Any realistic ribose-acid yield (0.3-1 per ribose x 35 mM consumed) breaks that at 4.5 and
    5.5 even with CBX at 0.
  - So Meynier tests the total acid against the buffer description: Davies at I ~0.6-1 M, plus Na+-P2O7 complexation
    that is not modelled. It is not a test of CBX specifically.
- **Hou.** The 20 min point is on a cliff: the 5 mM ARP is the only buffer.
  - ARP pKa 7.07 / 7.6 / 8.18 gives 6.35 / 4.48 / 3.88 (measured 5.05).
  - The trajectory-wide miss is the SHAPE: the model's acid stops when the ARP is gone, while the measured pH keeps
    falling 4.30 -> 3.55 at 60-120 min.
  - Consistent reading of both controls: the acid is allocated to the wrong source. The amino-acid carboxyl via
    AMINE_FATE carries it all; the sugar-fragment yield is fitted to about 0.
- **Do not just flip AMINE_FATE.** Carrying NH4+ passes Meynier but puts Zhou's endpoints at 5.0-5.9 against
  3.2-3.4 measured (AMINE_FATE_BASIS probe).
- **Verdict.** "Fails both": SUPPORTED, P 0.85. "CBX is the cause": WEAKENED, P 0.35.
- **Action.** Declared limitation and band now. Score unbuffered pots on measured trajectories where they exist. After
  A3: one arm that moves acid to sugar fragmentation, with a lagged NH4+ fate, scored jointly on
  Hou + Zhou + Martins + Meynier.

## Claim 2 (cluster 3): fed-H2S rows and the headspace

- **Kaw arithmetic.** Sander, B = 2100 K: 0.403 at 25 C and 2.17 at 145 C, checked by hand.
  - Real H2S volatility at 145 C (Fernandez-Prini 2003 / IAPWS kH ~185 MPa): Kaw ~1.0-1.1, close to the 100 C clamp
    (1.19).
  - So the dissolved fraction is about 0.24-0.26 and the bias is about 0.6 dex. The 0.84 dex end comes from
    van 't Hoff over-extrapolation.
- **Vessel.** Vg/Vl = 3 follows from "autoclave (200 mL)" with a 50 mL fill. That assumes 200 mL is the internal volume.
- **Speciation.** pH 5.34 in situ, alpha0 ~0.90: correct.
- **Do the fed rows ignore the partition?** Yes, and so does the WHOLE sulfur lane: `sulfur.py` has no H2S
  headspace term. The fed rows (`hofmann1998_t3/t4`) are in_model_fit and simulable.
  - So the fitted k is self-consistent for a Hofmann-type vessel, and biased only on transfer.
  - Cysteine-generated H2S partitions too: about 34 % dissolved in Hofmann's own 100 mL pots. A fed-rows-only fix
    would be inconsistent.
- **Unquantified.** H2S lost while preparing and charging the aqueous solution.
- **Verdict.** SUPPORTED, P 0.75 (bias 0.5-0.85 dex vs a dissolved-sulfide k).
- **Action.** A declared transfer band now. A global vessel operator for every H2S source after A3.

## Claim 3 (cluster 2): the hub

- **(a) Not separable: SUPPORTED, P 0.85, but it is structural.**
  - Quasi-steady, log L3/L1 = (lk1 - lk3 - a7) + (c1 - c3 - b)(pH - 7). Only those combinations are identifiable.
  - T8 shows it: A / B / N1 / N2 are all jointly infeasible (worst 0.33-0.44 dex) and all feasible per pot.
  - The negative controls N1 and N2 did NOT fail. The prereg's own rule therefore makes the test INCONCLUSIVE, not
    confirmatory.
  - The lk1 band (Martins, 0.49-49 /min) excludes Hou's own 1-DX bound (0.034-0.29 /min).
- **(b) Beck sets a band of 0.2-0.45 for pentose: WEAKENED, P 0.35.**
  - The source is a hexose, a propylamine ARP, two pH points, n = 1, GLC areas.
  - Buffer molarity is unprinted, the pH was not re-measured, and 0.2-1 M OPD (a base) was present for 10 h.
  - The dossier itself gives P(transfer) 0.3, and Martins' lower edge is also a hexose.
  - Counter-data: Zhou 2024 +-Ala has a measured dpH of 0.13-0.41 and a printed DR of about 0, against a predicted
    -0.03 to -0.25. In T9 the label spread persists (ratio 1.10).
- **(c) Drift is required: WEAKENED, P 0.5.**
  - Arm C fails only Hou per pot: 0.287 vs a tolerance of 0.25, with 6 parameters on 6 points. Every arm is jointly
    infeasible anyway.
  - Drift is measured, so the engine needs it as an operator. T8 does not show the ratios need it.
- **(d) No trap evidence survives: WEAKENED, P 0.3.**
  - Cystine is not a clean non-thiol: it is a 10x-supersaturated suspension, and in-situ reduction or exchange is
    flagged in the Hou dossier.
  - The clean non-thiol (Ala, 2 pots) gives DR ~1, and pH alone is rejected for the late cystine points.
  - So (iii) is diagnostic of a sulfur amino acid, not of a thiol. The cysteine DRs themselves vary from 0.45 to 1.8.
- **Action: nothing in the engine.** Declare the pentose split slope b in [0, 0.45] and the dpo_share as a band. The
  trap stays an uncredited arm, and Hou's cystine pot a weight-0.5 hold-out.

## Claim 4 (cluster 4): MFT loss and the observation operators

- **(a) Grosch: SUPPORTED, P 0.85.**
  - Verbatim, p. 189: "the addition of the deuterated internal standards immediately after the boiling period would
    compensate alterations caused by this reaction".
  - Methods: quantified "in the cooled reaction mixture (Method I) and in the SDE extract (Method II)".
  - Caveats:
    - The sentence is about the extract's storage oxidation.
    - Method II's standard goes into the extract, so SDE recovery is not compensated.
    - lambda assumes formation is equal under SDE and reflux. SDE strips H2S (Table VI, x3) and the precursor, so
      lambda is model-dependent and is probably a lower bound.
- **(b) The Nestle bound is NOT a measurement, P 0.15 that it is measurement-grade.**
  - It is derived from 4 single-vial, unreplicated, digitised points with co-heated standards of unprinted amount,
    through an assumed exhaustion profile.
- **(c) The law is FLEXIBLE, not identifiable, P 0.1 that it is right and pinned.**
  - Parameters: k_E, K_c, k_ox and a reversal rate, plus unmeasured E(t) and O2(t).
  - Data: 2 Grosch lambdas (one lab, thiamine pots, about 1 SD apart: the K_c CI includes no effect), 1 derived
    bound, and an O2 regression with n = 4 and p = 0.125.
  - The L0 rejection compares channels that differ (a thiamine-electrophile pot vs a co-heated patent). That violates
    the lesson "a second lab must measure the same channel".
  - P(pot-dependent >= 3x at matched T) ~0.5.
- **Operators.**
  - O1 / O2 / O3 are sound in kind, P 0.8. O1 caveat: a post-heat standard can exchange with the disulfide pool.
  - O4's rho 0.6 comes from one beef endpoint during SDE's own boil.
- **Action.** k_mft_decay as a declared band (Cys-conditional at most). The co-heated-standard tracer is a per-row
  operator after A3; until then the US 10,076,126 shape and level stay unscored. No law in the engine.
