# M1b pre-registration: the pot's charge-balance pH (part i) and the pentose-side pH laws (part ii) (v3 DRAFT)

**Status: v3.6 DRAFT; part (i) freezes after its power check.** No build, screen, fit or score of any rival exists. Part (i)
freezes when its screen and power check (secs. 3.5 - 3.6) have run and been answered; after any screen only the
demotion of a contrast to "undecidable on disk" is allowed, never a change to a band, a unit or a rule (sec. 8). Part
(ii) is a pre-declared conditional (sec. 4.0): its screen and power check run first, and its fits run only if the split
contrast survives.

*Plan step M1b ([STATUS_AND_PLAN](../../docs/guides/STATUS_AND_PLAN.md) sec. 5, W2, row M1b; sec. 7 items 1 - 9). Case:
[M1b_sugar_step_ph_sulfur_lane.md](../../docs/reference/cases/M1b_sugar_step_ph_sulfur_lane.md) (items 1 - 5 as drafted
2026-11-04, before the rescope; not edited here). Conflicts: CF-084, CF-060
([conflicts_register.yml](../../docs/reference/conflicts_register.yml)). Structural template:
[kinetic_core_m1_supply_loss_combinations_prereg.md](kinetic_core_m1_supply_loss_combinations_prereg.md) (v3.3, frozen),
its [results review](kinetic_core_m1_results_review.md), and the M1 + M2 draft (`kinetic_core_m1m2_supply_routes_refit_prereg.md`,
branch `m2-refit-prereg`, 9bc55883, unmerged) as amended by its adversarial review.*

- v1 2026-11-05 (cdb94a29); v2 the same day after an independent adversarial review (verdict SPLIT: part (i) freezes
  after fixes, part (ii) cut to what is decidable on disk) and the coordinator's decisions on it (9e649d0f); v3 after the
  re-review of v2 (MORE FIXES, then FREEZE AFTER POWER CHECK; sec. 19). No engine was run to write any version. "(arithmetic)" marks numbers computed here from printed or keyed values.
- Branch `m1b-prereg` from origin/main at cc5cbb55. `STATUS_AND_PLAN.md`, the registers and the case files are not edited;
  what this design needs from them is sec. 17. Every computation runs through `./scripts/docker_maillard.sh run`, inline
  pH coupling, shipped calibration B60.

## 0. Seen by the author

**The expectations in sec. 14 are not blind.** Read while writing:

- **Engine pH numbers in records.** CF-084 (cluster 3, 2026-11-01, **on B51 as shipped then, not on B60**): Meynier
  cysteine pots drift -0.79 / -0.65 / -0.39 (measured <= 0.2); Hou 2026 cooled 6.35 at 20 min -> 4.26 at 120 min against
  measured 5.05 -> 3.55, residuals pred - meas **+1.30 / +0.24 / +0.25 / +0.42 / +0.58 / +0.71 (cooled: the engine is too
  basic at every time)**; Zhou 2023 ARP alone 3.19 / 4.55 / 7.09 against 3.99 / 4.13 / 4.99; ARP pKa 7.07 / 7.6 / 8.18 RMS
  0.69 / 0.45 / 0.84 on Hou. These are re-printed at B60, cooled, before the screen (V6, sec. 3.4).
- **Cao 2017, corrected scale (coordinator, 2026-11-05).** v1 quoted "3.86 -> 3.89" against the measured 4.67 -> 4.15;
  that compared an in-situ value at 90 C with cooled readings. **On the cooled scale the engine starts at 5.500 against
  Cao's 5.50; under C1 it reads ~0.2 low at 30 min (4.476 vs 4.679), all of it CF-084's explanation 1 (7.9 mM of cysteine
  carboxyl carried with its amine declared destroyed).** The keying is correct; v1's audit item is withdrawn. Cao points
  toward B1. EX-013's Cao pH hold-out (RMS 0.080, entry arm with trunk solutes, after the forked-pool correction).
- W1.1's Meynier items (NF and FUR pH-flat to 14 digits in glycine + ribose); EX-033 at R3En (s042): Zhai 2021 Fig. 5
  pH-segment slopes of log10(3-DX / 1-DX), measured -0.21 / -0.38 / -0.80 and -0.14 / -0.38 / -0.70, predicted -0.02 /
  +0.40 / +0.09 and +0.01 / +0.18 / +0.08; EX-045 / EX-046; B61 Amendment 1; M1's carried optima; B60's fit report
  (`ph_drift`: acid yield 7.50e-5, on its lower bound; ARP pKa 7.068).
- **Measured values, no engine number:** Mikami 2017 Figs. 3A - C, 6B; Nakamura 2020 Figs. 2, 3, 7; Zhai 2021 Figs. 4b,
  5 and its LOD paragraph (PDF p. 10650); Hou 2026 Fig. 3 (ARP, pH, 3-DX, 1-DX, glycine); Davidek 2008 Fig. 1; Davidek
  2006a Fig. 3 (dossier); Davidek 2002's pH-stat acids (case); Meynier 1995's methods and results (PDF); Cao 2017 Fig. 4;
  Kang 2026's final pH; Zhou 2024 / 2025b's AX-ARP pH trajectories (CF-060); Martins 2003 Part II Table 3, Martins 2003b,
  Martins 2005b Table 3.
- The M1 + M2 draft, all of it. The adversarial review of v1 (its findings as relayed by the coordinator).
- **v3.4:** the dossiers `bork2022`, `bork2023`, `bork2022b`, `gupta2011`, `blank2006` (branch
  `lit-2026-11-05-norfuraneol-floors`, read in its worktree), including their engine arithmetic (the shipped NF sink and
  B61's term at pH 5 / 130 C).
- **v3.5:** the dossiers `sohn1995` and `carta1996` (branch `lit-2026-11-05-m1c-bands`, read in its worktree), and the
  Mizutani 1960 comparison inside the sohn1995 dossier.
- Not read: any engine number on Mikami's or Nakamura's lysine pots, Davidek 2008's or 2006a's pots, Hou's deoxyosones
  or glycine, Kang's final pH, Zhou 2024 / 2025b's pH, or any rival of this file.

## 1. Question, and the two separable parts

### 1.1 The claim

1. **(i) The pot's pH.** The charge balance (`ph_state`, B51) predicts in-situ and cooled pH from the network's
   titratable centres. It fails CF-084's controls in opposite directions: **too acid** in Meynier's buffered cysteine
   pots, **too basic** in Hou's unbuffered ARP pot at every time. Part (i) builds CF-084's two explanations (the
   released amine survives; the acid lags), the acid-supply rival the signed residuals point at, and their combination,
   and decides whether any balance may set the in-situ pH every pH law reads.
2. **(ii) The pentose-side pH laws.** In glycine (or lysine) + pentose pots the shipped engine's NF and furfural are
   pH-free and its 3- / 1-deoxypentosone ratio has the wrong sign against pH (EX-033, at R3En). v2 keeps only what the
   data on disk can decide (sec. 4.2): **the enolisation split carried to pentoses (P2) and a pH-dependent pentose +
   amine entry (P3)**. Base-catalysed NF decay (E16, P1) is declared **undecidable on disk**; the base-favoured
   3-deoxypentosone removal (P4) is screened, never fitted. **Not built:** `ph.trunk_terms` alone; E2's hexose part.

### 1.2 Separability, order, and what does not wait

- **Part (i)** has its own rivals, units, screen, power check, decision rule, stop rule, hour cap and ship path
  (sec. 3). It needs nothing from part (ii).
- **Part (i) does not gate the next thiol-supply case (M1c).** M1c runs at the measured pH trajectory where one exists
  (the measured-trajectory operator, MTO, sec. 3.1), whatever part (i)'s verdict. A balance part (i) validates later
  replaces the forced trajectory in those pots from then on.
- **Part (ii)** runs after part (i)'s verdict, on the pH of sec. 4.1, and is re-scored, no refit, under the other part
  (i) rivals (stacking guard, sec. 4.8). Part (i)'s controls are re-checked under part (ii)'s chosen rival (reverse
  stacking, sec. 9).

## 2. The evidence, by kind (case item 2; additions and corrections only)

The case's sec. 2 tables stand. Additions, each checked against the dossier, ledger, register or PDF named:

| rank | kind | source (lab) | what it shows | bears on | status |
|---|---|---|---|---|---|
| 2 | fed intermediates | Mikami 2017 Fig. 6B (Ochanomizu) | fed NF 1.72 mM, 0.2 M phosphate, **pH 8 only**, 100 C: k 0.032 /min | E16 at one pH | B61's anchor; a band source, never decisive |
| 2 | fed intermediates | Nakamura 2020 Fig. 3 | fed NF, pH 8 only | E16 at one pH | a data copy of Mikami's design (split v2 note "data copy (CF-019) of calibration input mikami2017") |
| 2 | fed intermediates | **Bork 2022 Fig. 1D (TU Berlin, Kanzler; v3.4)** | NF 20 mM **alone**, unbuffered water set to pH 5.0 +- 0.1 at room temperature (never re-measured: it drifts, probably down), sealed ampoules (headspace, degassing, chelator not stated), **130 C**: 100 / 85.9 / 75.1 / 59.2 / 42.2 / 20.0 % at 0 / 30 / 60 / 120 / 180 / 300 min; first order, **k 4.8 - 5.1e-3 /min** (triplicate; no LOD) | E16: a fed-NF point at a second pH, from a second laboratory; **and the shipped pH-free NF sink's size at acid pH** | dossier on branch `lit-2026-11-05-norfuraneol-floors` (`bork2022_extraction.md`; not yet on main); **a band source, not decisive** (sec. 4.2) |
| 2 | fed intermediates | Bork 2023 (same lab) | the same nominal pot, one endpoint: 48 +- 4 % left at 120 min, k ~ 6.1e-3 /min (two-point) | as above; with Bork 2022 the lab's range is **4.4 - 6.1e-3 /min**: one laboratory, the papers differ by more than either's replicate SD | `bork2023_extraction.md` (same branch) |
| -- | not used for E16 | Gupta 2011 (whey gel, pH 6.1 / 8.3, 105 C), Blank 2006 (xylose + Gly / Leu, pH 6, 70 C) | NF **formation**; no fed-NF arm ("NF decay not separable", gupta2011 dossier l. 36; Blank: endpoint only) | not E16's decay law | same branch; recorded |
| 2 | fed intermediates | Nakamura 2020 Fig. 7 | fed NF pH 8, 3 h: decay "strongly repressed" by 40 mM DTPA and by N2 bubbling (chromatograms) | **against** a pure hydroxide law: the decay is metal- and O2-dependent (rival P1') | qualitative; recorded |
| 2 | fed intermediates | Davidek 2006a Fig. 3 (Nestle) | synthetic 1-DG, **no amine**, 0.2 M phosphate, 120 C, 60 min: acetic acid 0.20 / 0.56 / 0.80 - 0.84 mol / mol at pH 5 / 6 / 6.5 - 8 (flat above 6.5); 1-DG gone within 5 min at pH >= 6 while the acid keeps rising | acid from deoxyosone cleavage, its yield and pH law (part (i) B3) | dossier [M] / [F]; **un-keyed** (adjudication 2026-10-30: nominal 1-DG concentration, overlapping pH 6.5 - 8 markers): a band source. **The printed yields are lower bounds**: the 1-DG concentration is nominal (deprotection and resin losses unprinted, so the true charge is at most the nominal), and at pH 5 the yield is still rising at 60 min. The cleavage is **downstream of 1-DG**: 1-DG is gone within 5 min at pH >= 6 while the acid rises 0.28 -> 0.56 at pH 6 over 60 min (dossier "use" item 3). **Same laboratory (Nestle) and the same SIDA method ([2H3]acetic acid) as the decisive Davidek 2008 units**: those units are not laboratory- or method-independent of B3's band. (The "pH 4 / 6 / 8" grid is Table 1's glucose + glycine pots; Fig. 3 is pH 5 - 8, 25 C labels.) Hexose 1-DG in ~0.16 M phosphate: a transfer to pentose 1-DX in other media |
| 5 | mass balance / measured product yield | **Sohn & Ho 1995 (ho_network, Rutgers; one of E10's DEV labs; v3.5)** | cysteine **alone** (no sugar), 0.05 M, water, **initial pH 8 (NaOH, unbuffered)**, 2 h, 0.3 L stainless-steel reactor (air presumed), NH3 by gas-sensing electrode: **~14 / ~76 / 82 % of the alpha-amino N released as NH3 at 110 / 150 / 180 C** (Fig. 2; only 82 % printed, the others figure reads); NH3 is a LOWER bound on the amine released and surviving (headspace partition, re-addition to carbonyls) | **FOR B1** (phi_Cys): the cysteine amine leaves largely as ammonia, a base, not destroyed. Bears on CF-084 explanation 1 and on Cao's ~0.2 cooled-pH offset (sec. 0) | dossier on branch `lit-2026-11-05-m1c-bands` (`sohn1995_extraction.md`; not yet on main). **Out of every part (i) objective, so a band source** (sec. 3.1); **not decisive** (not a CF-084 control). Caveats: one lab; unbuffered pH 8, not the controls' pH 4.5 - 6.5; steel reactor; cysteine only, so no sugar carbonyls to consume the NH3 (the very consumption `AMINE_FATE_BASIS` invokes) |
| 5 | mass balance | Mizutani 1960 (Hokkaido) | free cysteine 0.1 M, water, 120 C, open, CO2 sweep (acidic): NH3 0.75 % of N at 3 h, 1.34 % at 6 h | **conflict with Sohn of ~30 - 60x in the deamination RATE** (Sohn's 110 C first-order equivalent 1.3e-3 /min is ~30x Mizutani's ~4e-5 at 120 C; sohn1995 dossier sec. 5): base-favoured or vessel-catalysed cysteine deamination, separated by neither paper. It bears on the rate (unchanged in B1), not on the survival fraction phi_Cys; it does say Sohn's pH 8 yields need not transfer to acid pots | E10's TEST contrast F04#26 (not read for an engine number here); recorded |
| 6 | one factor in one study | Davidek 2002 (Nestle) | DFG 0.1 M, 90 C, **pH held** by a pH-stat: acid 3 vs 79 mmol/L at pH 5 vs 8 | acid supply rises with pH | evidence for B3's pH law; a **guard** (hexose, 90 C, hot electrode); not a band |
| 4 | measured intermediates | Hou 2026 Fig. 3A, 3H (Jiangnan) | GX-ARP 5 -> 3.20 / 2.18 / 1.27 / 1.10 / 0.75 / 0.68 mM (0 - 120 min); glycine released per ARP consumed 0.40 / 0.32 / 0.29 / 0.27 / 0.25 | the ARP's amino acid is released (B2's f_ret); the conditioning series for C-Hou | dossier [F] / [Z] |
| 4 | measured intermediates | Zhou 2024 Fig. 3b, Zhou 2025b Fig. 3d (Jiangnan) | the same nominal unbuffered AX-ARP 20 mM pot, initial pH 8, 120 C: pH 5.13 / 4.32 / 4.19 / 4.08 vs 6.04 / 5.05 / 4.66 / 4.50 (30 - 120 min); 3-DX 3.7 - 5.3 vs 0.37 - 0.52 mM; released Ala x2.2 apart | **the laboratory's own reproducibility of a pH trajectory: 0.4 - 0.9 pH units** between two papers of one design; the 3-DX scale ~10x across five AX-ARP pots | CF-060; reported checks (sec. 3.9), never controls; Zhou 2025b's dx pots are in part (ii)'s objective |
| 6 | one factor in one study | Martins 2005b Table 3 (Wageningen) | glucose + glycine condensation k1: pD 0.50 +- 0.14 dec / pH (95 % CI), 100 C, pH 4.8 - 7.5 | band of P3's entry slope | regressed on the **initial room-temperature label pH** of drifting pots (Martins 2005b sec. 2.1, "initial pH values ... measured at room temperature"): carried as a transfer (sec. 4.3) |
| 7 | fitted (engine record) | EX-003 (sub-model 1, arm C) | slopes on the glucose-consuming steps (the condensation among them): in-sample on its bound; **rejected on the frozen hold-out** (low-pH drifting runs 10 - 200x low in acid) | **against** a condensation pH law carried into the sulfur lane (P3) | engine record |
| 4 | measured intermediates | Lau 2003 (whey gel, 116 C), Pandit 2006 (mashed potato, 121 C) | flat NF plateaus (B61 audit item 5: 6 - 15 % loss at pH 6 - 6.5, 39 % at pH 7 predicted by B61 over Pandit's plateau) | **against** a strong hydroxide NF decay unless the matrices sit at pH <= 6.3 - 6.6 (unreported) | B61 audit; recorded |
| 4 | measured intermediates | Cao 2017 Fig. 4c (BTBU) | cooled pH, blank 5.49 -> 4.67 / 4.56 / 4.49 / 4.33 / 4.26 / 4.15 (0 - 5 h); x24 glycine arm 5.42 -> 4.43 ... 4.03; read +-0.03 | **a descriptive cross-check, not a control**: TTCA (ring N pKa 6.24) and CAMA titrate, so Cao's acid budget depends on M1c's supply choice. It points toward B1 (sec. 0) | development; seen |
| 4 | measured intermediates | Zhou 2023 Fig. 2 (Jiangnan) | Ala-xylose ARP 20 mM, water, 120 C, 60 min; final pH from labels 6 / 7 / 8: **ARP alone 3.99 / 4.13 / 4.99** (out of every objective); **ARP + Cys 3.22 / 3.42 / 5.07** (the three B60 fit rows `zhou_final_pH_from_pH6/7/8`, `ph_state` sec. "WHAT IS FITTED HERE"); read +-0.06 | CF-084's within-lab control | the ARP-alone arm is decisive (C-Zhou-ARP); the ARP + Cys arm a guard (C-Zhou-Cys) |

**Corrections to the plan row and the case (also sec. 17):**
1. "Fed norfuraneol loss vs pH" did not exist on disk through v3.3: both fed series were at pH 8 and one a same-lab
   repeat. v3.4: Bork 2022 / 2023 add NF alone at pH 5, 130 C, from a second laboratory; sec. 4.2 reassesses E16.
2. EX-033's wrong sign was measured at R3En (s042), not at B60.
3. The deoxyosone pH data are one laboratory group (Zhai 2021, Hou 2026: `ho_network`), semi-quantitative (CF-062), and
   Zhai 2021 Fig. 5 shaped P4 (EX-045 / 046).
4. Davidek 2008 is EX-013's fit data; no band or start here comes from it.
5. `ACID_YIELD_PER_SINK_EVENT` and `ARP_AMINE_PKA` are the balance's two fitted constants, **calibrated on exactly the
   three Zhou ARP + Cys rows** (`ph_state`, "WHAT IS FITTED HERE"); at B60 the yield sits on its lower bound (7.50e-5).
   They are integrated over declared bands in every part (i) rival (sec. 3.1), never held at their in-sample values.
6. **Every Martins pH slope** (2003 Part II, 2003b, 2005b) is regressed on the initial label pH of pots that drift (0.2
   - 0.9 units, Martins 2003 Part I sec. 3.1.6; up to ~1.4 in Martins 2005b's high-pH arm). v3.3 re-expresses each on
   the time-averaged pH (sec. 4.3.1), with the label-pH form as required agreement.

## 3. Part (i): the charge-balance pH

### 3.1 Rivals, coordinates, bands

All on the shipped structure (R0), `ph.b51` on, inline. The balance's two in-sample constants (rules 8 - 10):
- **The ARP secondary-ammonium pKa**, in every rival, uniform on **[7.07, 8.18]** (B60's value to Mossine 1994's
  measured DFG analogue, 8.18 +- 0.01, `ph_state` l. ~427: the incumbent and the only measured Amadori ammonium).
- **The lumped acid yield** of the deoxyosone sink (`r_osone_decay_*`, `r_arp_*` -> `ACID` x yield): in **B0, B1 and
  B2** uniform on **[0, 1]** (its meaning; `ACID_YIELD_BOUNDS`), and as required agreement **log-uniform on [1e-4, 1]**
  (B60's 7.5e-5 sits on its lower bound, so a uniform draw almost never visits the incumbent's region). In **B3 and
  B13 it is 0**: there the acid of the deoxyosone sinks is the measured acetic branch (below), not a lump, and keeping
  both would count the same acid twice. This resolves v2's contradiction (v2 said "every rival" and "set to 0").
- **B0-PINNED** (B60's shipped constants: yield 7.50e-5, pKa 7.068) is scored beside every rival: a favoured rival must
  also beat B0-pinned, as required agreement (sec. 3.7).

| rival | what changes | new coordinates | band (source) |
|---|---|---|---|
| **B0** | nothing (amine of released amino nitrogen destroyed at release, `AMINE_FATE_BASIS`; the ARP's amino-acid carboxyl carried as CBX at the enolisation) | -- | -- |
| **B1 amine survives** | the centres the ledger destroys on `AMINE_FATE_BASIS` keep a fraction as ammonium (pKa 9.25, the ledger's group), **separately** for cysteine thermolysis and its Strecker routes (**phi_Cys**) and for the Amadori / CAMA enolisations and `r_tdp_cys_trap` (**phi_ARP**) | phi_Cys, phi_ARP | each **[0, 1]** (fractions), primary. **v3.5: phi_Cys has a sourced form, required agreement:** **[0.76, 1]** from Sohn & Ho 1995 (NH3 >= 76 % of the cysteine N at 150 C, 2 h; since at most all the cysteine degraded, the surviving share of the released amine is at least the NH3 yield; at 110 C the 14 % bounds nothing because the degraded fraction is unknown). A transfer from a sugar-free, unbuffered pH 8 steel-reactor pot to the controls' sugar pots at pH 4.5 - 6.5, where released NH3 can be consumed by carbonyls: so it is a required-agreement form, not the primary band, and a B1 / B13 verdict counts only if it holds under both. phi_ARP: no source (stays [0, 1]) |
| **B2 acid lag** | the Amadori / CAMA enolisation returns its amino acid with fraction **f_ret** (as its zwitterion, no net centre) and destroys the rest as B0; the carboxyl equivalent B0 carries at ARP loss goes to a latent pool LCX (no charge) released as CBX first-order at **k_lag** | f_ret; log10 k_lag (120 C); Ea_lag | f_ret [0, 1] (Hou's 0.25 - 0.40 is unit C-Hou-Gly, not a band); k_lag from 0.0097 /min (Zhai 2023 JAFC TDP post-peak decline, 120 C) upward: Davidek 2006a's 1-DG lifetime (< ~2 min at pH >= 6, 120 C) gives k >= 0.35 /min, **a lower bound on the fast end, not an upper end**; the band's numerical upper end is 10 /min (lifetime 0.1 min, B0's timing to within the first sample), so a verdict conditional on it does not ship (rule 10); Ea the carbohydrate class distribution (B55), k(120 C) held |
| **B3 acid supply** | acid from deoxyosone cleavage at a measured yield, through a latent step (sec. 3.1.1): a share y(pH) of DPO's loss goes to a latent uncharged 2,4-diulose pool **LAA**, which releases **AA** (acetic acid, one carboxyl) first-order at k_iso; the lumped `ACID` yield 0 | y5, y6, y_hi (shares at in-situ breakpoints, sec. 3.1.1); log10 k_iso | **one-sided bands from the printed lower bounds** (Davidek 2006a Fig. 3, sec. 2): y_hi in [0.80, 1], y6 in [0.56, 1], y5 in [0.20, 0.56] (the upper end at pH 5 declared by monotonicity: every printed yield is non-decreasing in pH, so y5 <= the pH 6 value), with the order y5 <= y6 <= y_hi enforced per draw; each lower end widened by the **declared transfer term** delta_T (sec. 3.1.1); k_iso log-uniform on **[0.025, 0.24] /min at 120 C** (arithmetic: first-order approach fitted to Fig. 3's 5-min / 60-min yields, 0.03 / 0.20 at pH 5, 0.28 / 0.56 at pH 6, ~0.56 / ~0.82 at pH >= 6.5, giving k ~ 0.025, 0.14, 0.24), Ea 76.0 kJ/mol (`ea_dpo_aa_kj_mol`, Martins 2005's hexose analogue `r_odg_aa`); Davidek 2008 in no band |
| **B13** | B1 and B3 together | phi_Cys, phi_ARP, y, k_iso | as above; **lumped ACID yield 0** (as B3) |
| **MTO** (not a balance; the default) | the measured cooled trajectory converted to reaction temperature by the balance's own T-shift at the measured point (CF-084's inversion remedy); where no pH is measured, the **label pH converted to reaction temperature with the same shift** (not held at its 25 C value), with B0's drift as a sensitivity | -- | -- |

- **B3 makes acid through the `AA` species**, which today only `r_dpo_aa` (off) produces (`sulfur.py` l. 959; verified,
  no other step emits AA). So B3 also gives part (ii)'s acetic-acid units a source (sec. 4.2).

#### 3.1.1 B3's step, exactly

- **Code as it stands (verified 2026-11-05).** `r_dpo_aa` is DPO -> AA + FRAG_C 3 (`sulfur.py` l. 959); its rate is
  computed **outside** the entry's `if _ent.enabled` block (l. 1488 - 1489: `k_at(k_dpo_aa_per_min_100c_ph6,
  ea_dpo_aa_kj_mol, T) x aa_ph_factor(pH)`), and the step is switched off only when `k_dpo_aa_per_min_100c_ph6 <= 0`
  (l. ~997). `CENTRE_LEDGER` gives it one carboxyl (l. 1112). It is an **added** first-order DPO sink, not a branch.
- **The branch rule (declared): hold DPO's total loss, rescale the other branches.** At each right-hand-side call, with
  K_DPO the sum of DPO's first-order and pseudo-first-order loss constants in the rival's network at that state,
  the acetic branch gets y(pH) x K_DPO and every other DPO sink is multiplied by (1 - y(pH)). **DPO's lifetime is then
  unchanged in every rival**; what changes is where its carbon goes: every DPO-derived product (NF, DDP -> MFT, the
  C2 / C3 fragments) falls by (1 - y). That is the chemistry Davidek's yield implies (most 1-deoxyosone ends as acetic
  acid at pH >= 6.5) and it is a large side effect on MFT, recorded before the run (sec. 9). B3 therefore does not use
  the shipped `k_dpo_aa` / `aa_ph_factor` form; it uses the branch form with the band above. (The alternative, an added
  sink, shortens DPO's lifetime differently in each rival.) **Both forms are required agreement (v3.1):** the
  **added-sink form** puts the acetic branch on top of DPO's other sinks, at k_aa = y / (1 - y) x K_DPO,other at each
  state, so the share is still y but DPO's lifetime shortens by (1 - y); a B3 / B13 verdict counts only if it is the same
  under both forms.
- **The latent step.** The branch feeds LAA (uncharged; the 1-deoxy-2,4-diulose of Davidek's Fig. 7), LAA -> AA at
  k_iso. So acid appears after the deoxyosone is consumed, as Fig. 3 shows.
- **Breakpoints converted to in-situ pH.** Davidek's pH 5 / 6 / 6.5 are 25 C labels of ~0.16 M phosphate; each is
  converted to its in-situ value at 120 C by the balance's T-shift for that buffer at that label (printed by V5), and
  y(pH) is interpolated linearly in **in-situ** pH between them, flat outside.
- **Transfer term (declared): hexose 1-DG -> pentose 1-DX, ~0.16 M phosphate -> other media (unbuffered, pyrophosphate).**
  No hexose / pentose acetic ratio exists outside the decisive set (Davidek 2008 is decisive), so the transfer is a
  declared widening of each lower end by **delta_T = 0.30 dex** (a factor 2 downward). It widens downward only: the upper
  side is already open to 1 (the printed values are lower bounds). **Being unsourced, it must not decide the verdict
  (v3.1): required agreement across delta_T in {0, 0.30 dex}; a B3 / B13 verdict counts only if it holds under both.**
- **B12 is not built** (v1 had it): sec. 3.2 shows both its parts move Hou the wrong way.
- The case's third explanation (the ARP pKa) is now inside every rival through the pKa band.

### 3.2 What each rival predicts, from the signed residuals

| control | measured vs B0 (CF-084, B51; re-printed at B60 by V6) | B1 | B2 | B3 | B13 |
|---|---|---|---|---|---|
| Meynier cysteine pots (0.2 M pyrophosphate, Cys 50 + Rib 35 mM, 140 C, 1 h) | B0 **too acid** (drift -0.79 / -0.65 / -0.39 vs <= 0.2): 37.6 mM cysteine carboxyl, amine destroyed | right way (phi_Cys) | none (the acid is cysteine's) | wrong way (more acid), small where little 1-DX forms in a buffered cysteine pot (UNVERIFIED: the screen prints the 1-DX flux) | phi_Cys can offset B3 |
| Hou GX-ARP 5 mM, water, 120 C | B0 **too basic at every time** (+1.30 / +0.24 / +0.25 / +0.42 / +0.58 / +0.71; 6.35 vs 5.05 at 20 min) | **wrong way** (less acid) | **wrong way** (delays acid: the early miss grows) | **right way** (more acid, rising with pH while the pot is still near 7) | right way if phi_ARP is small |
| Hou level, 0 - 20 min | measured falls 1.95 (7.0 -> 5.05), B0 0.65 (-> 6.35) | wrong way | wrong way | right way | right way |
| Zhou ARP alone, labels 6 -> 8 | B0's endpoints spread 3.90 against 1.00 measured | depends on phi_ARP | flattens the spread if acid arrives late at high label | adds acid most at high label (y rises with pH): shrinks the spread | both |

So **neither CF-084 rival alone can pass both controls**: B1 fixes Meynier and worsens Hou; B2 worsens Hou. The signed
residuals point at **B13** (amine survival for cysteine, acid supply for the Amadori pot). The plan row requires
CF-084's two rivals to be built and compared; they are, and the comparison is expected to reject them alone (sec. 14).

### 3.3 Units, value basis and error terms

pH is read after cooling in every source; every rival is compared on its cooled (25 C) pH. **V5 (below) prints, for every
pH comparison, the scale (cooled 25 C, label 25 C, or in situ at T) and the temperature of both sides.**

**The pH error model (MAJOR 1).**
- **Measurement:** the printed replicate SD where one is readable; else max(read error, meter repeatability r), with
  **r = 0.05 pH (decided; sensitivity 0.02)**. A figure-reading error is not meter precision, so it never stands alone.
- **Model discrepancy, from H0's omega (rule 3), by exact per-species re-solves (v3).** **omega = 0.419 dex** (ln
  equivalent 0.965; EX-059: the worse of H0's two sets governs, H0 prereg sec. 3.6 / D2 and amendment A1 item 3; the
  warm set gives 0.419, the reset set 0.390; `kinetic_core_thiol_lolo_nested_tau.md` l. 23 and l. 80; the shipped H0
  application and the plan use 0.419) is a fractional error on what the network makes. At each draw and sampling time,
  the balance is re-solved twice with every **produced** titratable species (gross amounts formed, not the net: CBX,
  ACID, AA, FA, LCX-released CBX, surviving ammonium, returned amino acid) scaled **independently** by 10^(+omega) =
  2.62x or 10^(-omega) = 0.381x (arithmetic), the signs chosen so that one solve has every acid up and every base down
  and the other the reverse (pH is monotone in each species' amount, so these two corners bound every combination of
  independent +-omega scalings). **delta-pH = half the spread of the two re-solved cooled pH values.** Two balance
  solves per draw, no re-integration. Required-agreement sensitivities omega 0.5x / 2x (0.210 / 0.838 dex); **0.390
  (H0's reset set) a named, printed sensitivity only**, never required agreement (it lies inside the 0.5x - 2x range).
  This replaces v2's buffer-capacity linearisation (which used net acid and a local slope).

| unit | data | value basis | measurement noise | discrepancy | tau | floor | role |
|---|---|---|---|---|---|---|---|
| **C-Mey** (Cys pots, labels 4.5 - 6.5) | Meynier 1995: "no change of more than 0.2 pH unit occurred after heating" (all four amino acids, all five pH) | **censored**: likelihood P(\|d_pred + e\| <= 0.2), e ~ N(0, 2 r^2 + delta-pH^2) (two meter readings, before and after) | r | delta-pH | none | none | **control** (decisive); development; CF-084 values seen. **Replicate-draw rule** for the tie band: the bound's edge redrawn as 0.2 + N(0, 2 r^2) per replicate |
| **C-Mey (Gly pots)** | same statement | as above | r | delta-pH | none | none | **guard for B3 and B13, decided now**: their drift must stay within 0.2 (the acetic branch acts on the amine-free DPO flux, sec. 3.1.1, so B3 can acidify a glycine pot). For B0, B1, B2 every glycine pot predicts ~0 drift (no cysteine, the entry off, no ARP): a tie, no information |
| **C-Hou** | Hou Fig. 3B, GX-ARP alone, pH at 20 - 120 min (6 points) | absolute cooled pH, **with the measured ARP(t) imposed as an input** (Fig. 3A, interpolated monotone between the printed points): ARP's loss rate is forced to -dARP_meas / dt, distributed over the rival's own ARP branches in their predicted proportions at each state; everything downstream free. This is the MTO's construction applied to ARP instead of pH. Acid equivalents per ARP consumed (inverted from the cooled pH by the rival's titration, CF-084's remedy) printed beside it | max(read 0.05, r) | delta-pH | none | none | **control** (decisive); development; CF-084 values seen. The unforced run is a required-agreement form |
| **C-Hou-ARP** | Hou Fig. 3A, ARP 3.20 / 2.18 / 1.27 / 1.10 / 0.75 / 0.68 mM at 20 - 120 min, in the **unforced** run | C / C0 (calibrated UPLC-Q-TOF; charge 5 mM) | read +-0.1 mM (dossier), as sd | omega per point | none (C / C0) | none printed; declared none | **reported unit**: says whether a rival's pH miss is the balance's or the ARP decay's; never a claim or guard (the ARP decay is not part (i)'s mechanism) |
| **C-Hou-Gly** | Hou Fig. 3H, glycine released / ARP consumed (5 ratios) | ratio of two calibrated species of one pot | propagated read errors (+-2 mg/L glycine, +-0.1 mM ARP) | **2 omega^2 = 0.351 dex^2** (a ratio of two predicted points; sd 0.593 dex) | none | none printed (Tables S1 - S2 not on disk): B0, B1, B3 predict no glycine (the amino acid is destroyed): **structural zero**, "cannot represent", never summed | claim unit for B2 only |
| **C-Zhou-ARP** | Zhou 2023 Fig. 2, ARP alone, final pH at labels 6 / 7 / 8 | absolute cooled pH | max(read 0.06, r) | delta-pH | none | none | **control** (decisive); out of every objective; CF-084 values seen |
| **C-Zhou-Cys** | Zhou 2023 Fig. 2, ARP + Cys (3.22 / 3.42 / 5.07) | absolute cooled pH | as above | delta-pH | none | none | **guard only**: the three B60 fit rows; the balance's two constants were calibrated on them |

### 3.4 Verification before any score

- V1 every new key off reproduces B60's cost on inline (M1's V6 reference).
- V2 B2 with f_ret 0 and k_lag at its fast end reproduces B0 on every control to 1e-6 pH; V3 B1 at phi 0 reproduces B0;
  B3 at y 0 with the lumped yield at B60's value reproduces B0.
- V4 the strong-ion difference of every pot is constant at every node.
- **V7 (ARP forcing):** with forcing on, the rival's ARP reproduces Hou's interpolated ARP(t) to 1e-6 relative at every
  node; with forcing off, the run is bit-identical to the unforced rival. The forced loss is never negative (the
  interpolant is monotone); B3's branch rule (sec. 3.1.1) keeps DPO's total loss equal to the rival's without the branch
  to 1e-9 relative at every right-hand-side call.
- **V5 (scale):** every printed pH comparison carries, for both sides, its scale (cooled 25 C / label 25 C / in situ at
  T) and its temperature; a comparison of two different scales is a build error.
- **V6:** CF-084's engine numbers (Meynier, Hou, Zhou both arms) and Cao's re-printed at B60, cooled, inline, before the
  screen; the screen's expectations (sec. 3.2) are checked against them and any sign that differs is written into the
  record before any rival is run.
- V-density (sec. 5.1).

### 3.5 Mechanism screen (rule 13; >= 600 draws)

At B60, inline, scoring rtol 1e-8, **>= 600 Sobol draws per rival** over its integrated coordinates (including the two
balance constants); where a rival's feasible region is small (fewer than 10 draws pass), the fraction of the band that
passes and its bounding box are reported (coverage), and the draw count doubles once.

| sig | measured | B0 (expected, sec. 3.2) |
|---|---|---|
| sigma-Hou-level | pH falls 1.95 units in 0 - 20 min (7.0 -> 5.05) | 0.65 |
| sigma-Hou-shape | pH still falls between 60 and 120 min while ARP < 1.3 mM (4.30 -> 3.55) | CF-084: acid prompt with ARP loss |
| sigma-Mey | Cys pots: \|drift\| <= 0.2 at every label | fails at three labels |
| sigma-Zhou | ARP alone: endpoint rise from label 6 to 8 **smaller than B0's 3.90** (toward the measured 1.00): a sign | 3.90 (equal by definition) |

**Kill rule (one):** a rival is killed only if no draw reproduces all four; a would-be kill is confirmed at one fit (start
0, sec. 3.8). Passing is screening, never evidence. Cao's blank and the checks of sec. 3.9 are printed per draw as
descriptions. Artifact `kinetic_core_m1b_part_i_screen.{md,json}`.

### 3.6 Power check, then freeze (MAJOR 6)

Sec. 7's design on part (i)'s pairs (B1, B2, B3, B13 against B0; B13 against B3), from the screen's draws. A pair whose
power does not exceed its false-credit rate is demoted to "undecidable on disk" before the freeze. **Then part (i)
freezes.** After the screen, only such demotions are allowed (sec. 8).

### 3.7 Decision rule (part i)

Densities, tie band and adequacy gate: sec. 5.
1. **Comparison.** A rival is **favoured over B0** iff C-Mey (Cys) and C-Hou both favour it, C-Zhou-ARP does not favour
   B0, and the guard C-Zhou-Cys does not favour B0 after Holm (sec. 5.5); for B3 and B13 also the guard C-Mey (Gly) does
   not favour B0. C-Hou-Gly can add to a B2 verdict, never replace a control. **Required agreement (N2):** the same holds
   against **B0-pinned** (B60's shipped constants), and with the lumped acid yield of B0 / B1 / B2 drawn **log-uniform on
   [1e-4, 1]** instead of uniform on [0, 1].
2. **Validation for a class** iff favoured AND adequate (sec. 5.4) on every control of the class:

   | class | controls | covers by declared transfer |
   |---|---|---|
   | buffered amino acid + pentose (0.2 M phosphate / pyrophosphate), with or without cysteine | C-Mey (Cys) | Zhang 2025b, Davidek 2008, Mikami / Nakamura; Hofmann's 0.5 M phosphate pots (higher capacity: transfer declared) |
   | unbuffered Amadori | C-Hou, C-Zhou-ARP | Hou's and Zhou's design class |

   **Unbuffered TTCA pots (Zhai 2021 Fig. 4b / 5, Kang 2026) are covered by no control**: a validated unbuffered balance
   is applied to them only as a sensitivity beside the MTO until a TTCA-pot control exists.
3. Required agreement for every "favours" and "adequate": omega 0.5x / 2x; r 0.02; the unforced C-Hou; B0-pinned and
   the log-uniform yield (item 1); for B3 / B13 both delta_T values {0, 0.30 dex} and both branch forms (hold-total and
   added-sink, sec. 3.1.1); for B1 / B13 both phi_Cys bands ([0, 1] and Sohn's [0.76, 1], v3.5); the supply structure (R0 and M1's C1 carried optimum) for the
   pots carrying cysteine only (Meynier Cys, Zhou ARP + Cys; C1 changes nothing in a cysteine-free pot).

### 3.8 Fits (part i)

- **Primary: no refit**, each rival at B60's vector with its coordinates integrated (the lesson "run the change at the
  incumbent's constants with no refit").
- **Stacking guard:** one refit per favoured rival on M1's R0 objective (B60's 52 rows + the `k_mft_decay` prior), M1's
  optimiser settings and comparison rule (two starts, 1.92 nats, decisive-prediction agreement); starts M1's carried R0
  optimum with new coordinates at band centres / B8's perturbation / B60's vector at the B0-reproducing corner. Free:
  M1's R0 set (16 + `k_mft_decay`) plus the rival's coordinates; those the Laplace marks unidentified stay integrated.
  A favoured verdict counts only if it holds at the refit's carried optimum.

### 3.9 Reported checks (development; never used to choose)

Cooled pH under every rival, with V5's scale labels, on: **Cao 2017** (both arms, measured; points toward B1: its ~0.2
cooled-pH offset under C1 is exactly the "amine destroyed" bookkeeping of 7.9 mM cysteine carboxyl, and Sohn 1995 measures
that amine leaving as NH3); **Zhou 2024
and Zhou 2025b's AX-ARP pots** (measured; their 0.4 - 0.9 pH disagreement, CF-060, printed as the laboratory's own
reproducibility of a trajectory); Zhai 2021 Fig. 4b's eight pots (no pH measured); Zhang 2025b's blank (development, seen
under the engine this session; no pH measured); Kang 2026's three TTCA pots (measured final pH 5.1 / 4.9 / 4.7); Zhai 2023
Food Chem's "7.0 to 4.9" (text, `simulable: false`); Davidek 2002's held-pH acid amounts (a guard on B3's pH law, printed:
the rival's acid at pH 5 vs 8 for a fed DFG pot is not computable in the sulfur lane, so only its sign on a pentose ARP
analogue is printed). Each under R0 and C1.

### 3.10 Stop rule, outcomes, ship path, hour cap

**Stop rule.** After the screen, the power check and the stacking refits (no further rival, band or control): **if no rival
is favoured over B0 AND adequate on C-Mey (Cys), C-Hou and C-Zhou-ARP, the verdict is "balance not identified; the
measured-trajectory operator stays", and part (i) ends.** CF-084 stays open with the record; a new balance rival is a new
case.

The outcomes are exclusive, decided in this order:

| order | condition | verdict | ships |
|---|---|---|---|
| 1 | a rival validated for both classes | "balance validated" | that balance (kind 2 change to `ph_state` / `CENTRE_LEDGER`; `AMINE_FATE_BASIS` rewritten in the same PR) after the side effects, the results review and **owner sign-off at ship time** |
| 2 | else, a rival validated for exactly one class | "balance validated for <class>" | that balance for that class only; the other class keeps the MTO; same conditions |
| 3 | otherwise | the stop rule: "balance not identified; the MTO stays" | nothing |

Row 2 cannot co-occur with the stop rule: the stop rule fires only when row 2's condition fails.

**What follows.** Every pot with a measured pH trajectory runs under the MTO in M1c and part (ii) whatever the verdict. A
validated balance replaces the MTO for its class later, and the converted-label default in pots without a measured pH.
**Shipping part (i) moves the shipped pH everywhere in the sulfur lane** (every `base`, `acid`, `thiolate`, `hydroxide`
and H2S-speciation factor reads it): **owner sign-off is required at ship time.**

**If B3 or B13 ships, part (i) also moves shipped sulfur-lane chemistry, not only pH.** The acetic branch takes a share
y(pH) of DPO's carbon, so MFT and every other DPO-derived product (NF, DDP, the C2 / C3 fragments and what they feed)
change by up to a factor (1 - y): most at in-situ pH >= 6.5, where y >= 0.8 (before delta_T). Therefore, for a B3 / B13
ship: the side-effect record (sec. 9; plan sec. 7 item 7) **must include every thiol row** (B60's 52 rows by group, the
declared thiol hold-outs) **and the full held-out contrast benchmark**, before and after, under both branch forms; each
thiol loss names its candidate compensating error (sec. 9); and the owner sign-off is required **for this reason too**:
shipped MFT / FFT predictions move.

**Hour cap.** The cap is on **author work: ~14 h** (build, V-checks, screen, power check, answers, scoring, review
answers). **Compute is on top of it**: ~6 - 21 h wall (typical ~9.5 h) at <= 8 cores (sec. 13; v2's 7 - 22 h assumed every pot in every
draw, v3 counts the ~12 control pots per draw). At the author-work cap, whatever stage is reached, the stop rule's
verdict applies unless a rival is already validated.

## 4. Part (ii): the pentose-side pH laws, cut to what is decidable

### 4.0 A pre-declared conditional (v3)

1. **First, only the screen (sec. 4.4) and the power check (sec. 7)**, ~2 h wall (sec. 13).
2. **If the split contrast (P2) is demoted to "undecidable on disk" and only the entry contrast (P3) survives, part (ii)
   runs no fits.** P2 and P3 go to the joint M1c + M2 refit as fixed-band components (sec. 10), with the reason recorded:
   the entry contrast would only re-confirm EX-013 on the same Davidek 2008 data; and B3's own yield law already supplies
   0.17 of the 0.24 dex the split would need on U-Dav-Y, against a ratio discrepancy of 0.55 dex (the coordinator's
   arithmetic, 2026-11-05, from the re-review; not recomputed here: the screen prints both numbers, and if they differ
   materially the rule still applies, with the printed numbers as its reason).
3. **If the split survives** (with or without the entry), part (ii) freezes and runs as secs. 4.5 - 4.8.
4. If both are demoted, part (ii) ends with "undecidable on disk"; P2 and P3 are deferred to the joint refit.

### 4.1 Which pH part (ii) reads

The MTO wherever a trajectory is measured (Hou); a validated balance for its class; otherwise the label converted to
reaction temperature (MTO default). Zhai 2021's unbuffered TTCA pots are read **two ways, both required**: the converted
label, and the label drifting linearly to Kang's class endpoint (5.1 / 4.9 / 4.7 at 120 C, the class's only measured
endpoints), both in the power check (sec. 7).

### 4.2 What is decidable, settled before any run

- **The acetic-acid source, from the code (verified).** `AA` is produced only by `r_dpo_aa` (`sulfur.py` l. 959), off in
  B60. Its rate is computed outside the entry's block (l. 1488 - 1489), and DPO has an amine-free source,
  `r_pent_thermal_dpo` (`k_pent_thermal`, log10 -3.44 in B60's report). **So the step runs without the entry, but only
  at the amine-free DPO flux**: in P0 and P2 a glycine + xylose pot makes acetic acid only from thermal 2,3-enolisation;
  in P3 and P3+P2 also from PAMA's enolisation. Part (ii) carries B3's acetic branch (sec. 3.1.1: branch rule, latent
  step, Davidek 2006a band, never Davidek 2008) in every rival alike, so U-Dav-AA and U-Dav-Y are representable in each,
  with the amine-free-only flux as P0's / P2's honest prediction.
- **AA carries no charge in part (ii)** (N5): unless the validated part (i) balance is B3 or B13, the acetic branch is a
  product-only step whose AA is excluded from the charge balance (a build switch on `CENTRE_LEDGER`'s r_dpo_aa entry), so
  carrying it for the acetic units does not change part (ii)'s pH. If B3 or B13 shipped, AA is charged as validated.
- **E16 reassessed with Bork 2022 / 2023 (v3.4): still "undecidable on disk", and why.** Two fed-NF points now exist:
  Mikami, in-situ pH 8.08, 100 C, 0.2 M phosphate, aerobic (Nakamura's same-lab Fig. 7: the loss is suppressed by DTPA
  and by N2), k 0.032 /min; Bork, initial pH 5.0, 130 C, unbuffered, ampoule headspace unknown, k 4.4 - 6.1e-3 /min (one
  lab, two papers). The three laws, with the barrier unmeasured for NF alone and therefore integrated over the
  carbohydrate class distribution (B55's `ea_class_prior`, n = 20, range 30 - 237 kJ/mol, median 103.7) as required
  agreement: (a) **pH-flat** at Mikami's value projects to k(130 C) >= 0.066 /min at any positive barrier (Ea 30:
  x2.05, arithmetic), 13x or more above Bork; (b) **hydroxide only** (B61, n = 1) gives 3.2e-4 /min at pH 5 / 130 C at
  the median barrier (16x below Bork; arithmetic as the bork2022 dossier) and reaches Bork only at Ea ~218 kJ/mol
  (arithmetic: ln(4.9e-3 / 2.65e-5) x R / (1/373.15 - 1/403.15)), inside the class range; (c) a **partial order** (or
  flat + hydroxide) has n = 0.37 - 1.07 across the class band (arithmetic; the bork dossier's 0.55 - 0.76 is the 81.5 -
  143.5 kJ/mol sub-range). So on the face of it the pH-flat aerobic law is rejected and n > 0. **It does not decide E16,
  for four reasons that do not go away with more draws:** (1) **O2**: Mikami's loss is the aerobic, metal-dependent
  channel; Bork's ampoules may be O2-poor once the dissolved O2 (<= ~2 % of the 400 umol lost) is spent, the headspace
  being unstated (the bork2022 dossier's O2 budget: an equal-volume air headspace would carry <= ~50 %), so the two points
  may measure different channels (P1') and the "pH" contrast may be an O2 contrast; (2) **the barrier**: the 30 C gap
  (100 vs 130 C) is bridged only by an unsourced-for-NF class band, which by itself spans n = 0.37 - 1.07, i.e. every
  law except pH-flat; (3) **Bork's pH drift**: unbuffered and never re-measured, NF breakdown making acids; a lower mean
  pH lowers n further (at a mean pH of 4.5, n = 0.31 - 0.92, arithmetic), and no in-situ pH at 130 C is printed;
  (4) **NF's self-reaction with its own fragments (v3.6, Bork's supplementary files, branch `lit-2026-11-05-bork-si`)**:
  MGO and glycolaldehyde, NF's own fragments, are strong NF sinks at pH 5 (a 20 mM partner adds roughly 0.012 - 0.16 /min
  over the first 30 min, Bork's binary pots), so part of Bork's k at 20 mM NF may be second order in NF and its
  fragments. That part scales with NF concentration and weakens any transfer to Mikami's 1.7 mM or to pot-level NF
  (micromolar to low millimolar); Bork ran one NF concentration, so it cannot be checked on disk. The supplementary
  files add no headspace, O2, final pH or temperature series; their NF-alone series is Bork 2022 Fig. 1D replotted, with
  replicate SDs giving k = 4.91 +- 0.23e-3 /min.
  Separating hydroxide-only from partial-order needs a fed-NF temperature series at one pH, or fed NF at two pH in one
  laboratory with O2 controlled; neither exists on disk. **Bork's role is therefore a BAND, not a unit**: it sets the
  size of the NF loss at acid pH for the joint refit (sec. 10), and it is decisive for nothing in M1b.
- **Independent of the pH law: the shipped pH-free NF sink is too slow at pH 5.** The shipped `r_nf_decay` (log10
  k(145 C) -2.182, Ea 143.5) gives 1.41e-3 /min at 130 C (arithmetic); Bork measures 4.4 - 6.1e-3: **3.1 - 4.3x the
  engine** (3.4x at Bork 2022's 4.9e-3), with B61's term adding only 3.2e-4. Evidence against the shipped constant at
  acid pH whatever law governs pH; recorded for the register (sec. 17) and handed to the joint refit (sec. 10), not acted
  on in M1b (E16 owns no change here).
- **E16 (P1, NF decay pH law), as declared through v3.3:** one fed pH (Mikami 8, Nakamura a copy); the pot NF shapes
  (Mikami Figs. 3A - B, Nakamura Fig. 2) are confounded by the entry level, lysine's reactivity with NF, and O2 / metals
  (Nakamura Fig. 7). The oxidative rival **P1'** (metal / O2-catalysed loss through the enolate) is named; no quantitative
  datum separates it from P1. Lau 2003 / Pandit 2006's flat plateaus are evidence against a strong hydroxide law. **Any
  later E16 statement is scoped to aerobic, unchelated pots.** P1 is screened (sec. 4.4) and its units are reported
  (sec. 4.6), never fitted or credited. **Precondition for any future E16 credit:** the entry's NF yield passes s6 (the
  EX-013 entry arm made 25 - 36 mol % NF against 4 - 16 % measured).
- **P4 (3-DX removal): screen and signature only**, no fits (its claim units are the data that shaped it).
- **What remains fitted and judged:** P0, P2 (split), P3 (entry), P3+P2.

### 4.3 Rivals, coordinates, bands, starts, the free parameter set

| rival | steps changed | new coordinates | band (source; never a scored unit) | starts 0 / 1 / 2 |
|---|---|---|---|---|
| **P0** | none, plus B3's acetic branch (all rivals; uncharged unless B3 / B13 shipped, sec. 4.2) | -- | -- | -- |
| **P2 split** | the amine-free pair (`r_pent_thermal_dpo` vs `r_pent_caramel`) and, with the entry on, PAMA's pair: odds(1-DX : 3-DX) x 10^(s (pH - 5)), total flux unchanged; the cysteine-catalysed pair keeps its `base` / `acid` factors | s (an **odds** slope, dec / pH, on the in-situ pH it is applied on) | **[-0.02, 0.62]** (v3.3; drift-corrected, sec. 4.3.1): the range of measured hexose odds slopes, each = slope(1-DG route) - slope(3-DG route), re-expressed on the time-averaged pH: Martins 2003 Part II 100 C 0.55 x 0.810 = 0.45, 120 C -0.03 x 0.790 = -0.02; Martins 2005b k7 - k4 0.25 x 1.22 - 1.47 = 0.31 - 0.37; Beck 1988 ~0.44 (no drift reported: kept on its label basis, flagged); Martins 2003b's standing 3-DG : 1-DG (12 at pH 5.5 vs 3.2 or 1.2 at 6.8; the same pots as Part II) 0.44 - 0.77 x 0.79 - 0.81 = 0.35 - 0.62 (arithmetic). **The label-pH band [-0.03, 0.77] is a required-agreement form** (the 5.5 arm's drift shape is text-only). Hexose -> pentose: transfers declared | 0.30 / perturbed / -0.02 |
| **P3 entry** | `PENTOSE_AMINE_ENTRY` on: `r_pent_amine` x 10^(e (pH - 6)); PAMA's enolisations at the shipped factors (4.14 / 9.13 x ARP's); `amine_return` 1 | log10 k_entry (100 C, pH 6); e | k_entry fitted on the in-sample Desclaux 2006 pot, bounds +-3 dex about the Desclaux-probe value (UNVERIFIED on disk; else Martins' hexose k_schiff, 96.8 kJ/mol, transferred); **e ~ N(0.67, 0.18^2)** (v3.3; drift-corrected, sec. 4.3.1): Martins 2005b's k1 pD 0.50 re-expressed on the time-averaged pH, x 1.22 - 1.47 = 0.61 - 0.73, centre 0.67; sd = sqrt(0.095^2 + 0.145^2 + 0.06^2) = 0.18 (arithmetic): the 95 % CI's sd 0.071 x 1.34, the declared transfer term for hexose -> pentose and glycine -> lysine (the sd of Martins 2003 Part II's six per-step slopes, themselves drift-corrected: 0.300 / 0.745 / 0.583 at 100 C, 0.545 / 0.521 / 0.600 at 120 C, sd 0.145; 0.18 on the label basis), and half the 0.61 - 0.73 spread from the drift filter. **The label-pH form e ~ N(0.50, 0.19^2) (v3.2) is a required-agreement form** (Martins 2005b's drift is read at +-0.15 pH). EX-013's 0.57 is not used | probe value / perturbed / -1 dex; e 0.67 / 0.49 / 0.85 |
| **P3+P2** | both | as both | as both | as both |
| P1 (screen only) | `r_nf_decay_oh`: k8 x 10^(n (pH_insitu - 8.082)) | n, Ea_oh | k8 0.032 /min +- 10 % at Mikami's in-situ pH 8.082 (100 C); n [0, 1]; Ea the class distribution | -- |
| P4 (screen only) | `r_tdp_x`: TDP -> FRAG_C (no MGO, no thiol), x 10^(m (pH_insitu - 5)) | k_x, m, Ea_x | k_x from off (-10) to the measured total TDP loss 0.0097 /min at 120 C (Zhai 2023 JAFC), **converted to pH 5 at that pot's in-situ pH over the decline window** (from the balance or MTO, printed), per draw's m; m [0, 1]; Ea the class distribution | -- |

**Free parameter set of every part (ii) fit:** M1's R0 set (B60's 16 identified coordinates and `k_mft_decay` under its
prior) + P2: s; + P3: log10 k_entry, e; P3+P2: all three. Held at B60: everything else, B3's acetic shares at their band
draws (integrated), the balance constants integrated as in part (i) unless part (i) shipped a balance. Starts never from
scored data.

#### 4.3.1 Martins' slopes re-expressed on the pH they are applied on (v3.3; as M1c and M3)

Every Martins slope is regressed on the **initial room-temperature label pH** of pots that drift while they heat. P2 and
P3 apply their slopes to the in-situ pH at each state, so a slope is converted by the ratio of the label gap to the
**time-averaged** gap (a rate integrates over the trajectory). The drift shape inside a run is not printed for every arm,
so time-averaging assumes a linear drift; the label-pH band is kept as a required-agreement form in both rivals.

| source | drift used, and its source line | label gap -> time-averaged gap | factor on the slope | applied to |
|---|---|---|---|---|
| Martins 2003 Part II (and Part I = Martins 2003b, the same pots), 0.1 M phosphate, DFG | initial 6.8 arm -0.19 after 180 min at 100 C and -0.21 after 120 min at 120 C; initial 5.5 arm -0.8 and -0.9 (text only, no time course): `martins2003_extraction.md` l. 103 and l. 407 - 410; `martins2003b_extraction.md` l. 96 and l. 704 - 712 | 1.30 -> 1.605 (100 C) / 1.645 (120 C) (arithmetic: 6.8 - 0.095 vs 5.5 - 0.40; 6.8 - 0.105 vs 5.5 - 0.45) | 0.810 / 0.790 | P2's Part II and 2003b components; P3's transfer term (the six per-step slopes) |
| Martins 2005b, 0.1 M phosphate, glucose + glycine, 100 C, five arms 4.8 - 7.5 | Fig. 4 table, read +-0.15 pH: final at 240 min ~4.65 / 5.25 / 5.55 / 5.9 / 6.1 (`martins2005b_extraction.md` l. 725 - 751); the paper keeps only points with a drop <= 1 unit (l. 116, l. 245 - 246), so the pH 7.5 arm's effective mean is "nearer ~7.0" (l. ~748). **Opposite asymmetry: here the high-pH arms drift most**, so the time-averaged gap is SMALLER than the label gap and the slope RISES | regression of the label on the trapezoid time-averaged pH over 0 - 240 min: 4.72 / 5.36 / 5.72 / 6.24 / 6.53 (or 7.0 for the filtered 7.5 arm) (arithmetic) | 1.47 (full window) / 1.22 (filtered 7.5 arm) | P3's e (k1); P2's Martins 2005b component (k7 - k4) |
| Beck 1988 | phosphate (concentration unprinted), pH not re-measured, 10 h reflux with OPD: **no drift reported** (`beck1988_extraction.md` l. 60) | -- | none (label basis, flagged "drift unknown") | P2 (0.44 kept) |
| B3's breakpoints | taken from **Davidek 2006a, not Martins**: no Martins correction applies (their conversion to in-situ pH is sec. 3.1.1). B3's Ea 76.0 kJ/mol (Martins 2005's hexose `r_odg_aa`) is not a pH quantity | -- | -- | -- |

**Old -> new bands:** P2's s [-0.03, 0.77] -> **[-0.02, 0.62]**; P3's e N(0.50, 0.19^2) -> **N(0.67, 0.18^2)**; the
transfer term inside P3's sd 0.18 -> 0.145. Both old forms are required agreement (sec. 4.8's stacking list includes
them), and both are fitted as separate engine forms only where a fit runs (sec. 4.0's conditional); in the screen and
the power check they are separate draw sets.

### 4.4 M1b.0: signature screen (rule 13; >= 600 draws)

At B60, inline, part (i)'s pH; >= 600 Sobol draws per form (P0, P1, P2, P3, P4, P3+P1, P3+P2, P3+P4) with sec. 3.5's
coverage rule.

| sig | data | measured | role |
|---|---|---|---|
| s1 | Zhai 2021 Fig. 5, both ladders: d log10(3-DX / 1-DX) / d pH per segment | negative in all six | claim sign for P2, P4 |
| s2 | Hou Fig. 3C - D: log10(3-DX / 1-DX) rises as the measured pH falls | positive in time | sign guard |
| s2' | Hou Fig. 3C: 3-DX falls 102.5 -> 65.5 mg/L over 60 - 120 min while pH falls 4.30 -> 3.55 | 3-DX removed at pH ~4 despite residual ARP | **sign guard on P4**: a base-favoured removal must still reproduce a net 3-DX fall as the pot acidifies |
| s3 | Mikami Fig. 3A: NF at 60 min has an interior maximum (pH 6.5) | interior maximum | reported (E16 undecidable) |
| s4 | Mikami Fig. 3B: NF falls 1 -> 2 h at pH 7, 8, rises at pH 6 | sign change | reported |
| s5 | Davidek 2008: acetic(pH 8) > acetic(pH 6); xylose consumed faster at pH 8 | positive | claim sign for P3; **a form that cannot represent it (no AA or no xylose consumption) drops from the contrast, it is not killed** |
| s6 | NF yield scale, Mikami / Nakamura 4 - 16 mol % (EX-013's entry arm made 25 - 36 mol %) | -- | **precondition for shipping P3 (N8) and for any E16 credit**: P3 (or P3+P2) ships only if its carried optimum's NF yield on Mikami's xylose + lysine pot at pH 6.5, 60 min, lies inside the measured 4 - 16 mol % range widened by its own predictive interval; printed at the screen, never a kill |

**Kill rule:** a fitted form (P2, P3, P3+P2) is killed only if no draw reproduces every claim sign (s1 for P2, s5 for P3,
both for P3+P2) and no guard (s2, s2') is violated by every draw; confirmed at one fit. P1 and P4 are screened for the
record: P4's s1 / s2' pass fractions are its whole result in M1b.

### 4.5 Objective: rows in and out

| block | rows | measurement model |
|---|---|---|
| B60's 52 rows | B60's | B60's |
| `k_mft_decay` prior | M1 sec. 2.6 | as M1 |
| pentose block | `deoxyosone_submodel.yml`: `desclaux2006_xyl_gly_ph6_100c` (P3's k_entry); `zhou2025b_axarp*_dx` (the pentose-ARP split at one drifting pH: P2's level); `zhai2023jafc_*_dx` at 120 C (not P4's: P4 is not fitted; kept so that P0 and the fitted rivals share one objective) | read error / printed SD + s_series (M1's protocol) |

**Out of every objective:** every unit of sec. 4.6; Mikami Fig. 6B and Nakamura Fig. 3; Zhai 2021 Figs. 4b, 5; Davidek
2006a (band only), 2008 (decisive); Meynier 1995 (item 7); Hou 2026; Zhou 2024 (CF-060's other pot; reported).

### 4.6 Units (rules 1 - 3)

omega **0.419 dex** (ln 0.965; sec. 3.3) per predicted point in every primary Sigma: a ratio carries 2 omega^2 =
0.351 dex^2, sd 0.593 dex (arithmetic; 0.552 at v3's 0.390); contrasts sharing a reference point the corresponding
covariance. Declared transfer of omega from thiols to these species; required-agreement sensitivities 0.5x / 2x (0.210 /
0.838 dex) and within-pot rho 0.5; 0.390 a named, printed sensitivity only. No unit carries tau (each is a ratio or a
within-study contrast).

| unit | data | value basis | measurement noise | floor (primary; sensitivity) | role |
|---|---|---|---|---|---|
| **U-Zhai-T**, **U-Zhai-TX** | Zhai 2021 Fig. 5a - b: log10(3-DX / 1-DX)(pH, t) minus pH 5.5's, t 60 - 120 min (12 contrasts per ladder) | contrast of a two-species ratio of one chromatogram (shared diacetyl IS) | printed SD (n = 3) with read error | 1.50 umol/L: LOD 0.306 mg/L for both quinoxalines (PDF p. 10650) / 204.2 g/mol (arithmetic); read as the pot; LOQ 5.20 umol/L | claim (P2); "consistent with" only (shaped P4) |
| **U-Hou-DX** | Hou Fig. 3C - D: log10(3-DX / 1-DX) at 20 - 120 min relative to 20 min | within-pot time contrast; pH and time confounded | read +-2 / +-0.7 mg/L as sd | transfer: Zhai 2021's quinoxaline LOD; LOQ | guard only |
| **U-Dav-AA** | Davidek 2008 Fig. 1: log10 AA(pH 8) / AA(pH 6) at 30 / 60 / 120 / 240 min | within-study ratio, SIDA | read error with a declared CV 10 % (Davidek 2006a / b's; 2008's own unstated) | none printed: primary none (floored count printed); sensitivity 1 mol % (the papers' reporting resolution) | claim (P3) |
| **U-Dav-X** | Davidek 2008: ln(Xyl / Xyl0) at pH 6 and 8, common times (pH 8 points on the axis excluded as keyed) | C / C0 (one assay, the charge known) | read error | transfer Zhai 2021's Xyl LOD 0.020 mmol/L; LOQ 0.070 | claim (P3) |
| **U-Dav-Y** | Davidek 2008: AA / (Xyl0 - Xyl), pH 6 vs 8 | branching ratio | as above | as above | claim (P2) |
| U-Mik-pH, U-Mik-t | Mikami Fig. 3A, 3B (Xyl + Lys, 0.2 M phosphate, 100 C) as v1 | within-panel ratios | read +-3 % of scale; cross-panel +-15 % | transfer: Zhai 2021's DAD LOD as molar, 2.7 umol/L; x10 | **reported** (E16 undecidable) |
| U-Nak | Nakamura Fig. 2, pH 6.5 and 8: per-pH time shapes **and the cross-pH ratio NF(8) / NF(6.5) at 0.5 / 1 / 3 / 6 h** (both panels are in mM: v1's reason for dropping it was false) | within-study ratios | the dossier's read error | as U-Mik | **reported** unless the power check shows it decides a contrast of sec. 4.7 (then a guard) |

Lysine (Mikami, Nakamura) is charged as the entry's amine, 2x charge as a sensitivity (reported units only now).

### 4.7 Contrasts

| contrast | pairs | claim units | guard units |
|---|---|---|---|
| split (E2, pentose) | P2 vs P0; P3+P2 vs P3 | U-Zhai-T, U-Zhai-TX, U-Dav-Y | U-Hou-DX |
| entry | P3 vs P0 | U-Dav-AA, U-Dav-X | U-Hou-DX |

The split ships only with U-Dav-Y favouring it (non-Zhai support); Zhai-only credit is recorded as "consistent with".
Removal and E16: no contrast (sec. 4.2).

### 4.8 Decision rule and outcomes

Sec. 5's rules; a contrast credits a iff its claim units' block favours a in every listed pair, no guard favours b, and the
Holm in-sample guard holds. **Stacking guard:** each credit must hold, no refit, under every part (i) rival's balance and
under the MTO, and under both Zhai pH readings; otherwise "conflict: the pH decides it". **Required agreement on the
bands (v3.3):** each credit must also hold with P2's and P3's label-pH bands ([-0.03, 0.77]; N(0.50, 0.19^2)) in place of
the drift-corrected ones (sec. 4.3.1).

| result | ships |
|---|---|
| split credited with U-Dav-Y | P2's s posterior (provisional until the joint refit) |
| entry credited | P3 (on in every glycine / lysine pot: a chemistry change; owner sign-off), **only if s6 holds** (sec. 4.4); else "credited, held: NF yield precondition failed" |
| Zhai-only split credit | nothing; recorded |
| none / ties / both inadequate | nothing; the case's "cannot decide" list stands with the power check's numbers |

## 5. Densities, tie band, adequacy, multiplicity (rules 3 - 6, 9, 11)

### 5.1 Predictive density, per rival, then differenced (rule 4)

For unit u and rival k, lppd_uk = log of the integral of the unit's likelihood over the rival's parameter distribution
(Laplace on identified coordinates; profiles on bound-active ones; **integrated coordinates drawn from their bands**).
**Each rival's lppd is computed separately and then differenced**: d_u = lppd_u(a) - lppd_u(b).

- **Estimator:** Pareto-smoothed importance sampling with a unit-tuned proposal (multivariate Student-t, 5 df, at the
  mode of likelihood x posterior found on the unit's pots, <= 300 pot evaluations; scale the inverse central Hessian;
  band coordinates uniform in their bands). Accepted iff **k-hat < 0.7 and ESS >= 100** and ESS meets the tie-band
  precision rule (SE_MC(lppd_uk) <= 0.15 sd_rep(d_u) / sqrt 2). S starts at 256, doubles to 4 096; then Gauss-Hermite
  quadrature on the coordinates carrying > 90 % of the linearised variance; else "density not resolved" (a tie).
- **Exact densities on every form** (v1's restriction to forms that can flip is withdrawn). **Forms that change only the
  likelihood** (noise scales, r, floors, the censoring edge, rho, delta_T's effect is NOT one of them: it changes a band)
  **reweight the same engine draws**, with k-hat and ESS re-checked per form; a form whose reweighting fails the check
  gets its own draws. omega's 0.5x / 2x need only the two balance re-solves per draw (no re-integration). **Engine-
  changing forms get new runs:** the supply structure (R0, C1; cysteine pots only), the forced vs unforced C-Hou, the
  yield's log-uniform measure and delta_T = 0 (new draws of the integrated coordinates; the pots are re-run), B3's
  added-sink branch form (B3 / B13; v3.1), Zhai 2021's
  two pH readings, the part (i) balances under which part (ii) is re-scored, the lysine charge (reported units).
- Plug-in and linearised densities printed as diagnostics; a unit whose linearised value differs from PSIS by more than
  its tie band is flagged. **V-density:** the estimator reproduces M1's known case (R0's RES_Cys, alternative floor, ~ -746
  by quadrature) within its MC SE before any score.

### 5.2 Tie band (rule 5)

sd_rep(d_u) from R = 4 000 redraws of y_u from its **measurement noise only** (never omega; for C-Mey the edge-redraw rule
of sec. 3.3), at fixed predictive distributions; **b_u = z_0.95 sd_rep(d_u)** (the plan's 90 % convention).

### 5.3 Required agreement

A unit favours a iff d_u > b_u under the primary density and keeps its sign under every required-agreement form of secs.
3.7, 4.1, 4.6 and 4.8. Otherwise a tie.

### 5.4 Adequacy gate (rule 6)

A unit may favour (or veto for) a only if a is adequate there (posterior-predictive p > 0.05 for the Mahalanobis
discrepancy under Sigma_u including omega / delta-pH, for C-Mey the censored likelihood's PPC) or the rivals' predictive
means differ by more than the unit's discrepancy scale (omega x sqrt(c) for log units; delta-pH for pH units). Otherwise
**"both inadequate"**, listed with both p values.

### 5.5 In-sample guard with Holm (rule 11)

Per row group of the objective: Delta c_g and its one-sided p from the rows' declared sigmas; **Holm's step-down at
family-wise 0.05 one-sided** across the pair's groups. A lost group must name a candidate compensating error (sec. 9) and
must not outrank the best claim unit won; otherwise "credited, held".

## 6. Compensating structure (rule 10)

Bound-active coordinates profiled at three points inside the bound (0.25 / 0.5 / 1.0 dex; 1 / 3 / 10 kJ/mol; fractions and
orders 0.05 / 0.1 / 0.2 inside); switched-off rates (-10) get one extra point at the other rival's value; each decisive
d_u recomputed per point. A verdict conditional on a **numerical** bound (k_lag's 10 /min end; any -10 floor; the +-3 dex
fit bounds) does not ship; measured-class and definitional bounds (fractions 0 / 1, orders 0 / 1) do. B60's acid yield is
bound-active at 7.5e-5: part (i) integrates it, and B0-pinned is printed beside B0.

## 7. Synthetic power check (rule 7)

For each pair of secs. 3.6 and 4.7: "truth" = the screen's best draw of each rival; R = 200 synthetic data sets per truth
from the full Sigma_u (measurement + omega / delta-pH); both rivals scored with the primary density using the screen's
draw ensembles reweighted as the parameter distribution (no refit: overstates power); the full rule applied. Reported:
power, false credit, P(both inadequate). **A pair whose power does not exceed its false-credit rate is "undecidable on
disk"** and is not run; power below 0.5 is labelled "low power" in every verdict sentence. Part (ii)'s check runs under
**both** Zhai 2021 pH readings and a pair counts as decidable only if it is decidable under both. Budget: no new engine run
beyond the screens' draws; ~0.3 h on 2 cores per part.

## 8. Order of operations

0. H0's omega (on disk); M1's carried optima (on disk).
1. Review of v2 (sec. 12 (a')), answers.
2. **Part (i):** build B1, B2, B3, the MTO and the scorer, off by default; V1 - V6 and V-density; screen; power check;
   **freeze part (i)**; no-refit scoring; refits of favoured rivals; scoring at the refit; side effects; results review;
   decision (or stop rule); EX entry; owner sign-off if it ships.
3. **Part (ii):** build P2's split law, B3's branch form and the uncharged-AA switch (P3, P1 exist; P4's step for the
   screen only), off by default; V-checks (each switch at its P0 corner reproduces P0 to relative 1e-9); screen; power
   check; **the conditional of sec. 4.0** (no fits if the split is demoted); otherwise freeze part (ii); fits
   (two at a time at 4 workers, memory checked, never paused, `failed_pot_runs` read, M1's A1 rule); central polish;
   Laplace; profiles; scoring after every fit has ended; side effects; reverse stacking (sec. 9); results review;
   decision.
4. **After any screen, only the demotion of a contrast to "undecidable on disk" is allowed, never a change to a band,
   unit, form or rule.** No decisive number is computed on an intermediate member.

## 9. Side effects and reverse stacking

Before (R0-shipped (b)) and after (each shipped or favoured rival), every row: the full held-out contrast benchmark
(coverage first), with Meynier 1995's NF and FUR ladders read here only; B60's 52 rows per group; route residuals by pH;
the checks of sec. 3.9; the Zhou 2023 pH 6 / 8 columns (declared hold-outs). **Reverse stacking:** part (i)'s controls
(C-Mey, C-Hou, C-Zhou-ARP) re-scored, no refit, under part (ii)'s chosen rival; a part (i) verdict that changes is reported
and part (ii) does not ship until it is answered.

| expected loss | why | candidate compensating error |
|---|---|---|
| `zhou_final_pH_from_pH6/7/8` (C-Zhou-Cys) | the balance's two constants were calibrated on them; B1 / B3 move them | the destroyed-amine declaration and the near-zero acid yield standing in for each other on three rows |
| every MFT row and `hofmann_ribose_NF_insitu` (B3, B13; and part (ii)'s carried acetic branch where it diverts DPO carbon) | the branch rule sends a share y(pH) of DPO to acetic acid: DPO-derived NF, DDP -> MFT and C2 / C3 fall by (1 - y), most at pH >= 6.5 | `k_dpo_nf`, `k_dpo_ddp`, `k_nf_mft` fitted with no acetic branch, i.e. under a DPO pool that keeps all its carbon |
| `zhou_pH7_*`, `kang_*` thiol rows | every pH factor reads the new pH | B60 constants fitted under B0's pH |
| `hofmann_ribose_NF_insitu`, `_FUR_insitu` (P2) | the amine-free split moves | `k_pent_caramel`, `k_tdp_fur` fitted under pH-free amine-free steps |
| Kang / Zhou furfural rows (P2, P3) | the pentose's fate moves | the single furfural barrier (E11, M2) |
| Kumazawa FFT retention | none expected | a loss is a build error until shown otherwise |

**Tags:** Martins 2005b's k1 pD spent as P3's band; Mikami Fig. 6B B61's anchor; Zhai 2021 Fig. 5 "shaped P4"; Davidek
2008 "EX-013 fit data"; Davidek 2006a Fig. 3 spent as B3's band (un-keyed); Cao 2017 pH "EX-013 hold-out".

## 10. Interaction with M2 and the joint M1c + M2 refit

- **The dependency.** The M1 + M2 draft's Q2 (C1F-NF) uses a pH-flat NF decay with Mikami's pH 8 rate as its ceiling at
  every pH. M1b owns E16's pH law and **declares it undecidable on disk** (sec. 4.2).
- **Bork 2022 / 2023 as the joint refit's band at acid pH (v3.4).** Any NF-decay rival of the joint refit (Q2, Q2-pH)
  is banded at pH <= ~5.5 by Bork's measured k(130 C) 4.4 - 6.1e-3 /min (one lab), anchored at 130 C with the barrier
  integrated over the carbohydrate class distribution, scoped **"at ~20 mM NF, unbuffered, initial pH 5, O2 not
  stated"**; the aerobic pH 8 ceiling (Mikami) is not carried to pH 5 pots. **Concentration-transfer caveat (v3.6):**
  the band's upper end may include second-order self-reaction of NF with its own fragments (sec. 4.2, blocker 4), which
  does not transfer to pot-level NF; a first-order sink banded on it is an upper bound at lower NF. Bork is then a band in the joint refit, never one of its decisive units.
  The 3.1 - 4.3x shortfall of the shipped pH-free sink at pH 5 (sec. 4.2) is the joint refit's to resolve.
- **What enters the joint refit.** P1 is **declared deferred**. Q2 stays as written, with its verdict sentence carrying
  "the ceiling is Mikami's aerobic, unchelated pH 8 value; its value at the pot's pH is unmeasured (M1b)". **Q2's
  pH-dependent version, if any later step builds it, never carries an aerobic k8 into sealed or degassed pots** (Hofmann's
  autoclave, Cerny pot D's vials): Nakamura Fig. 7 makes the aerobic value an upper bound there. Where the joint refit
  prints the deferred rival's effect on its decisive unit (Cerny pot D's fed-NF share), the pH factor is 10^(n (8.082 -
  pH_insitu,potD)) with **pot D's in-situ pH at 95 C from the balance (or MTO), printed by the build**: v1's "10^(n x
  3.08)" subtracted a label pH from an in-situ one and is withdrawn.
- **P4** (screen only) moves the deoxyosone pool upstream of NF; it is deferred to the joint refit with its effect on pot
  D's share printed there at the screen's passing draws (m 0 / 0.5 / 1, k_x at its converted ceiling), no refit.
- **P2 / P3**, if credited before the joint refit is frozen, enter it as fixed-band components (their M1b posteriors).
  If part (ii)'s conditional (sec. 4.0) stops before any fit, they enter as fixed-band components with their **prior
  bands** (s in [-0.02, 0.62]; e ~ N(0.67, 0.18^2); the label-pH forms [-0.03, 0.77] and N(0.50, 0.19^2) beside them,
  sec. 4.3.1) and the recorded reason; otherwise they are deferred with their effect
  on pot D's share printed the same way. B3's acetic branch, if part (i) ships B3 or B13, is part of the shipped balance
  the joint refit inherits; it moves pot D's DPO carbon by (1 - y) at pot D's in-situ pH, printed by the build.
- **M1c** does not wait for part (i) (sec. 1.2); part (ii) is re-scored under M1c's supply winner if it exists by then.

## 11. Differences from the case and the plan

1. Fed NF across pH does not exist; E16 is undecidable on disk; P1 deferred.
2. Part (i) is separable, capped by a stop rule and an hour cap, judged on CF-084's controls (Meynier, Hou, Zhou ARP
   alone; Zhou ARP + Cys a guard), does not gate M1c, ships only with owner sign-off. CF-084's "RMS <= 0.3" and "within
   0.2" are replaced by measured noise, delta-pH from omega and the adequacy gate (Meynier's 0.2 is its own bound).
3. **B3 (acid supply) is added** and **B12 dropped**: the signed residuals (Hou too basic everywhere) say both CF-084 rivals
   move Hou the wrong way.
4. Case R2 - R5 map to P2, P1, P3, P4; only P2 and P3 are fitted; R1 (`ph.trunk_terms`) not built.
5. Furfural has no admissible pH datum (Meynier item 7; Zhou 2023 pH 6 / 8 declared hold-outs).
6. P4 is not EX-046's structure (no MGO, no thiol) and is screened only.

## 12. Adversarial review

(a) Design review of v1: done (SPLIT; findings in sec. 19). (a') **A re-review of v2** by a fresh agent before part (i)'s
build: the strongest case that B3 or B13 is built to win (Davidek's hexose yields, the integrated balance constants), that
the ARP forcing hides an error, that B3's branch rule or one-sided band is built to win, that delta-pH's re-solve
corners are too wide or too narrow, and that the outcome order lets a partial
validation through. (b) Results reviews per part before its decision; not sec. 14. Every objection is answered in the
record or the part waits.

## 13. Compute (current engine, inline, <= 8 cores)

Basis (`profile_fit_evaluation_2026-11-05.md`): R0's objective 3.92 s serial, 1.67 s at W = 4; median pot 0.068 s, slowest
0.46 s; the machine delivers **~3.5 - 4 serial-equivalents** at any worker count, so wall = serial-equivalent / ~3.5;
M1's R0 fits took 146 - 280 evaluations in 309 - 718 s at W = 4. **Part (i) per draw runs only the control pots: ~12
pots** (5 Meynier Cys, 1 Hou forced + its unforced run, 6 Zhou), ~0.2 s each (profile: Zhou pots 0.16 - 0.25 s),
**~2.6 s per draw**; B3 / B13 add the 5 Meynier Gly guard pots (~3.6 s); **C1 is run only on the 8 cysteine pots** (~2 s).
The descriptive checks of sec. 3.9 run once per rival at its carried point, not per draw. Part (ii)'s screen ~25 pots
(~3.75 s per draw). Every job `nice`d; parallel fits only if both fit in RAM.

| step | serial-equivalent | wall |
|---|---|---|
| **part (i)** author work | -- | **capped at ~14 h** (sec. 3.10); everything below is compute on top of it |
| screen: 5 rivals x >= 600 draws x ~2.8 s (mean) | ~2.3 h | ~0.7 h |
| power check | ~0.3 h | ~0.3 h |
| unit-tuned mode searches: ~5 units x 5 rivals (plus C1 on the cysteine units) x <= 300 x ~0.25 s | ~0.8 h | ~0.25 h |
| PSIS, B0 / B1 / B2: 3 rivals x S 1 024 - 4 096 x (primary ~2.6 s + log-uniform yield ~2.6 s + C1 on cysteine pots ~2 s) | 6 - 25 h | 1.8 - 7 h |
| PSIS, B3 / B13 (v3.1): 2 rivals x S 1 024 - 4 096 x (primary ~3.6 s + delta_T = 0 ~3.6 s + **the added-sink branch form ~3.6 s** + C1 ~2 s) | 7.5 - 29 h | 2.1 - 8.4 h |
| PSIS, Sohn's phi_Cys form (v3.5): B1 / B13 x S 1 024 - 4 096 x the 8 cysteine pots (~2 s) | 1.1 - 4.5 h | 0.3 - 1.3 h |
| (likelihood-only forms reweighted; omega's scales by the two balance re-solves per draw) | | |
| stacking refits: 1 - 2 rivals x 3 starts x 150 - 300 evaluations x ~3.9 s | 0.5 - 2 h | 0.3 - 1 h |
| checks (sec. 3.9) and side effects (for a B3 / B13 ship: every thiol row and the full benchmark under both branch forms) | -- | 1 - 4 h (UNVERIFIED scorer time) |
| **part (i) compute** | **~18 - 64 h** | **~6 - 21 h wall (typical ~9.5 h), on top of the author work** |
| **part (ii) conditional first step** (sec. 4.0): screen 8 forms x >= 600 x ~3.75 s + power check under both Zhai readings | ~5.5 h | **~2 h** |
| **part (ii) fit path** (only if the split survives), below | | |
| **part (ii)** author work | -- | 8 - 12 h |
| screen: 8 forms x >= 600 x ~3.75 s | ~5 h | ~1.4 h |
| power check (two Zhai readings) | ~0.5 h | ~0.5 h |
| unit-tuned mode searches: ~8 units x 4 rivals x ~3 engine forms x <= 300 x ~1 s | ~8 - 12 h | ~2.5 - 3.5 h |
| fits: 4 forms x 3 starts x 150 - 300 x ~4.3 s | 2.2 - 4.3 h | 0.7 - 1.3 h |
| central polish + continuations | 0.7 - 1.3 h | 0.2 - 0.4 h |
| Laplace | ~0.3 h | ~0.1 h |
| profiles: 2 - 6 coordinates x 3 points x <= 600 x ~4.3 s | 2 - 13 h | 0.6 - 3.7 h |
| PSIS: 4 rivals x ~3 engine forms (Zhai 2 readings; part (i) balances) x S 1 024 - 4 096 x ~3.75 s | 13 - 51 h | 3.7 - 15 h |
| side effects | -- | 1 - 3 h |
| **part (ii) compute** | **~32 - 88 h** | **~10 - 29 h (typical ~16 h); about 2.5 - 3.5 calendar days with the author work** |

Mode searches across both parts: **~10 - 14 h serial** (the coordinator's ~14 h). The densities dominate; the exact
density on every form (no restriction) is the cost of rule 4, partly recovered by reweighting likelihood-only forms.

## 14. Expectations (written now; not blind)

| event | P |
|---|---|
| V6: B60's cooled Hou residuals keep the B51 sign (too basic at every time) | 0.85 |
| part (i) screen: B1 alone fails sigma-Hou-level; B2 alone fails it | 0.85 / 0.85 |
| part (i) screen: B3 passes; B13 passes | 0.45 / 0.5 |
| part (i): B13 or B3 favoured over B0 | 0.35 |
| part (i): validated for the buffered / unbuffered Amadori class | 0.25 / 0.15 |
| part (i): the stop rule fires | 0.6 |
| power check: a part (i) pair demoted | 0.35 |
| power check: the entry or split contrast undecidable under at least one Zhai reading | 0.55 |
| part (ii): split credited with U-Dav-Y / entry credited | 0.1 / 0.2 |
| part (ii): no credit (the case's "cannot decide" stands) | 0.6 |
| part (i) verdict by 11-12; part (ii) by 11-20 | 0.5 / 0.2 |

## 15. Owner-level decisions

1. **Shipping part (i)** moves the shipped pH everywhere in the sulfur lane: **owner sign-off at ship time** (recorded;
   nothing to decide now).
2. **P3 turns the pentose + amine entry on** in every glycine / lysine pot (a chemistry change): ship on the rule, or hold.
(Decided and closed: the meter repeatability, 0.05 with 0.02 as the sensitivity; delta_T, required agreement across
{0, 0.30 dex}; B3's branch rule, required agreement across the hold-total and added-sink forms. An unsourced choice
must not decide the verdict.)

## 16. Open questions

1. *(Answered in v3 from the code, sec. 4.2: `r_dpo_aa` runs without the entry, at the amine-free DPO flux.)*
2. Does `TRUNK_SOLUTES_IN_BALANCE` cover a returned alanine (B2)?
3. Is the Desclaux-probe entry value on disk as a number (P3's start / bounds)?
4. ARP forcing (C-Hou) imposes the measured loss on the rival's branches in their predicted proportions; a rival whose
   branch proportions are themselves pH-dependent (all of them) sees a forced loss at a pH it did not produce. Accepted
   as the MTO's analogue; the unforced run is required agreement.
5. Kang 2026's final pH may be Zhai 2023 Food Chem's (CF-005 covers the SI only).
6. Zhai 2021's LOD basis (pot vs derivatised solution) is unstated.
7. Davidek 2008's own replicate statement (10 % CV transferred).
8. Davidek 2006a Fig. 3's pH 6.5 / 7 / 8 markers are ambiguous; B3 draws y_hi over 0.80 - 0.84 regardless.
9. Should part (ii) wait for M1c's supply winner?

## 17. For the register session

**CF candidates (ids to be assigned):**
1. *Fed norfuraneol at one pH only* (Mikami 2017 Fig. 6B, Nakamura 2020 Fig. 3: one design, pH 8).
2. *Cao 2017's pH on the cooled scale* (coordinator's numbers): the engine starts at 5.500 against 5.50 and reads ~0.2 low
   at 30 min under C1 (4.476 vs 4.679), all of it CF-084 explanation 1 (7.9 mM cysteine carboxyl, amine destroyed): append
   to CF-084 as a cross-check pointing to B1. v1's "3.86 -> 3.89" compared an in-situ value with cooled readings: withdrawn.
3. *The balance's two constants are calibrated on the ARP + Cys arm of Zhou 2023 only*, and the acid yield sits on its
   lower bound at B60 (7.50e-5): register the bound-active calibration.
4. *`AMINE_FATE_BASIS` vs Hou 2026's glycine*: 0.25 - 0.40 glycine released per GX-ARP consumed contradicts "the amino acid
   destroyed" in `r_arp_dpo` / `r_arp_tdp`.
5. *Zhai 2021 Fig. 5 is both P4's source and its only test* (spent-data note).
6. *Nakamura 2020 Fig. 7*: NF decay is metal- and O2-dependent; B61 encodes specific-base hydrolysis.
6b. *Bork 2022 / 2023 vs the shipped NF sink (v3.4)*: NF alone at pH 5 (initial, unbuffered), 130 C, k 4.4 - 6.1e-3 /min
   (one lab) against the shipped pH-free `r_nf_decay` 1.41e-3 /min and B61's 3.2e-4: the shipped constant is 3.1 - 4.3x
   too slow at acid pH, independent of the pH law. (The bork2022 dossier proposes the CF text; the dossiers are on branch
   `lit-2026-11-05-norfuraneol-floors`. Naming hazard noted there: the PDF on disk `bork2022.pdf` belongs to the key
   `bork2022b`.)
6c. *Sohn & Ho 1995 vs `AMINE_FATE_BASIS` (v3.5)*: cysteine alone at initial pH 8 releases >= 76 % of its amine N as NH3 at
   150 C (82 % at 180 C): the cysteine amine is released as a base, not destroyed; append to CF-084 (explanation 1). And
   the Sohn vs Mizutani 1960 deamination-rate conflict (~30 - 60x; pH 8 sealed steel vs acidic open CO2 sweep) as its own
   conflict entry (dossier on `lit-2026-11-05-m1c-bands`).
6d. *Carta 1996: Hou 2026's cystine pot is ~7x over saturation, not ~10x.* L-cystine saturates at 0.69 - 0.71 mmol/L in
   water at 25 C, pH 3 - 8.6 (carta1996 dossier, same branch); Hou charges 5 mM, so 7x. The hou2026 dossier's "~0.11 g/L,
   ~0.5 mM, my value" (sec. 1, answer 1) should be corrected to Carta's measured 0.168 g/L, 0.70 mM.
7. *Zhou 2024 vs 2025b pH trajectories* already in CF-060 (0.4 - 0.9 pH apart): add "the laboratory's reproducibility of
   an unbuffered pH trajectory" as its consequence for any pH control from this group.

**Wording errors found:**
1. Plan row M1b, decisive data: "Fed norfuraneol loss vs pH": no pH contrast exists in fed NF.
2. Plan row M1b and the brief: "wrong sign ... (EX-033)" is an R3En result, not B60.
3. Case sec. 2.2: Mikami is B61's anchor, so it cannot be decisive for E16.
4. Case sec. 2.7: Martins 2005b's k1 slope (0.50 +- 0.14) is missing; every Martins slope is on initial label pH.
5. CF-084 remedy: round-number tolerances ("RMS <= 0.3", "within 0.2").
6. EX-059's conclusion prints omega(all) 0.42 (warm) and the H0 review both 0.4187 / 0.3901: the governing value is the
   worse set's, 0.419 (H0 prereg sec. 3.6 / D2, A1 item 3), as the plan and the shipped H0 application use; M1's and
   the M1 + M2 draft's 0.390 should be aligned (v1 - v3.1 of this file used 0.390).
7. Davidek 2006a's amine-free 1-DG yields are at pH 5 / 6 / 6.5 - 8 (Fig. 3), not pH 4 / 6 / 8 (that grid is Table 1's
   glucose + glycine pots).

## 18. Results

*(appended per part after its run; nothing above this line is edited after the first number of that part exists)*

### 18.1 Part (i), 2026-11-06: built, screened, power-checked; CLOSED FOR SHIPPING

- **Records** (branch `m1b-part-i-build`, PR #94): `kinetic_core_m1b_part_i_checks.{json,md}` (BIT, V1 - V7,
  V-solver, V-density: all pass), `kinetic_core_m1b_part_i_screen.{json,md}` with its draws
  `kinetic_core_m1b_part_i_screen_draws.json.gz` (0 failed pot runs; no would-be kill), and
  `kinetic_core_m1b_part_i_power.{json,md}` (power 0.000 and false credit 0.000 for all five pairs: all demoted).
  Part (i) was never frozen; no `--fit` or `--score` ran.
- **Verdict (the coordinator's decision, 2026-11-06, after an independent review): "balance not identified; the
  measured-trajectory operator stays".** Part (i) is closed for shipping. The reason is NOT "the data cannot decide":
  1. **The +-omega corner error term (sec. 3.3, N3) was a mis-transfer.** omega (0.419 dex) is H0's between-pot
     scatter of thiol CONCENTRATIONS. Scaling every produced titratable pool independently by 10^(+-omega), acids up
     and bases down together, breaks the stoichiometry the pools share: Meynier's cysteine pots would carry ~98 mM of
     cysteine carboxyl (37.6 mM x 2.62) from 50 mM of cysteine charged. The two corners are a joint worst-case bound,
     and they were used as a 1-sd scale. The result was delta-pH of 0.4 - 3.8 pH, which (the review's finding)
     flattens B1's and B13's predictive densities about 4x for structural reasons.
  2. **CF-084's controls are spent.** They nominated every rival (sec. 3.2 was written from their signed residuals),
     so under the plan's sec. 7 (item 2, as the coordinator's decision cites it) they can falsify a rival but never
     credit one.
- What remains is the one falsification-only amendment of sec. 20 (A1, v3.7). The v3.6 power record stands beside it.

## 19. Revision log

- v1 DRAFT, 2026-11-05 (cdb94a29): first version, with the coordinator's corrections of the same day (Cao a descriptive
  cross-check; part (i) separable, stop rule, does not gate M1c, owner sign-off; part (ii) feeds the joint M1c + M2 refit
  or is deferred).
- v2 DRAFT, 2026-11-05, after the independent adversarial review of v1 (SPLIT) and the coordinator's decisions, before
  any build:

| finding | change |
|---|---|
| Cao scale | v1 compared in-situ 90 C with cooled readings; corrected to the cooled numbers, the audit item removed, Cao points to B1 (secs. 0, 2, 3.5, 3.9, 17); V5 prints both sides' scale and temperature |
| CF-084 numbers from a B51 record | re-printed at B60, cooled, before the screen (V6) |
| MTO "label held" | the label converted to reaction temperature with the balance's shift (3.1) |
| P4's ceiling | converted at the pot's in-situ pH (4.3) |
| sec. 10's 10^(n x 3.08) | withdrawn: label vs in-situ; the build prints pot D's in-situ pH (10) |
| Martins' slopes | stated as regressed on initial label pH of drifting pots, carried as transfers (2, 4.3) |
| BLOCKER 1: sec. 3.2 not from signed residuals | rewritten: Hou too basic at every time; B1, B2 move it the wrong way; B12 dropped (3.2) |
| acid supply rival missing | B3 (Davidek 2006a's 1-DG yields through `AA` / `r_dpo_aa`; Davidek 2002 a guard; Davidek 2008 in no band) and B13 (3.1) |
| in-sample balance constants held | integrated over [0, 1] and [7.07, 8.18] in every rival; B0-pinned printed (3.1, 6) |
| B2 one coordinate; one phi_N | B2 = f_ret + lag; phi_Cys and phi_ARP separate (3.1) |
| level signature | sigma-Hou-level added (3.5) |
| comparison not conditioned on ARP | C-Hou time-warped to the measured ARP(t); acid per ARP consumed printed (3.3) |
| MAJOR 1: pH error term | measurement max(read, declared meter r); discrepancy from omega through buffer capacity, per draw; C-Mey censored with a replicate-draw rule (3.3) |
| C-Zhou mixed in-sample and out | split: C-Zhou-ARP decisive, C-Zhou-Cys guard; sigma-Zhou a sign (3.3, 3.5) |
| stop rule and outcome table overlapped | an ordered exclusive table; hour cap ~14 h (3.10) |
| MAJOR 6: freeze before power check | part (i) freezes after its power check; after any screen only demotions (3.6, 8) |
| MINORs | C-Hou-Gly 2 omega^2; Meynier's glycine pots tie; B2's 0.35 /min a lower bound; s5 drops, not kills; s6 an E16 precondition; >= 600 draws or coverage; Zhou 2024 / 2025b (CF-060) in sec. 2; reverse stacking (9) |
| part (ii) too large | representability settled from the code (AA only from `r_dpo_aa`: B3's step carried in every rival); E16 undecidable, P1' named, scoped to aerobic unchelated pots, P1 deferred, Q2-pH never carries an aerobic k8 into sealed pots, Lau / Pandit listed; U-Nak's cross-pH ratio restored (reported); P2 band re-derived as an odds slope [-0.03, 0.77]; P3 widened by a declared transfer term, EX-003 listed; P4 screen only with Hou's 3-DX sign guard; free parameter set named; both Zhai pH readings in the power check (4, 7) |
| compute | exact densities on every form; per rival then differenced; likelihood-only forms reweighted with k-hat re-checked; mode searches budgeted; totals recomputed (5.1, 13) |

- v3 DRAFT, 2026-11-05, after the re-review of v2 (MORE FIXES, then FREEZE AFTER POWER CHECK) and the coordinator's
  decisions, before any build:

| finding | change |
|---|---|
| N1: B3's band | one-sided bands from the printed lower bounds (nominal 1-DG; pH 5 still rising), y5's upper end by monotonicity; a declared transfer term delta_T 0.30 dex downward (sensitivity 0); a latent 2,4-diulose step with k_iso [0.025, 0.24] /min from Fig. 3's 5 / 60 min yields; the branch rule (hold DPO's total loss, rescale the others; DPO's lifetime unchanged in every rival); breakpoints converted to in-situ pH; Davidek 2006a declared same lab and SIDA method as the decisive Davidek 2008 (secs. 2, 3.1.1) |
| N2 | favoured rivals must also beat B0-pinned and hold with the yield log-uniform on [1e-4, 1]; B13's lumped yield 0; the "every rival" vs "set to 0" contradiction resolved (3.1, 3.7) |
| N3 | delta-pH by two exact balance re-solves per draw with every produced species scaled independently by 10^(+-omega), gross, half-spread (3.3) |
| N4 | C-Hou with Hou's measured ARP(t) forced as an input (V7); ARP(t) reported as C-Hou-ARP; the time warp withdrawn (3.3, 3.4) |
| N5 | in part (ii) AA is product-only and uncharged unless B3 / B13 ship (4.2) |
| N6 | verified in the code: `r_dpo_aa`'s rate outside the entry block (l. 1488 - 1489), an added sink switched off only by its rate; DPO's amine-free source `r_pent_thermal_dpo`: the acetic branch runs without the entry at the amine-free flux (3.1.1, 4.2; open question 1 closed) |
| N7 | C-Mey (Gly) a guard for B3 / B13, decided now (3.3, 3.7) |
| N8 | s6 a precondition for shipping P3 (4.4, 4.8) |
| compute | part (i) per draw on ~12 control pots, C1 only on the 8 cysteine pots: ~5 - 17 h wall on top of the ~14 h author-work cap (3.10, 13) |
| meter repeatability | 0.05, sensitivity 0.02: decided; removed from the owner items (3.3, 15) |
| part (ii) | a pre-declared conditional: screen + power check first (~2 h); no fits if the split is demoted, P2 / P3 to the joint refit as fixed-band components with the reason (4.0, 10, 13) |

- v3.1, 2026-11-05 (coordinator's decisions on v3's two new owner items; no other change): delta_T is required agreement
  across {0, 0.30 dex}; B3's branch rule is required agreement across the hold-total and added-sink forms (the extra
  engine form budgeted, sec. 13: part (i) compute now ~17 - 59 h serial, ~6 - 20 h wall); both removed from the owner
  items (an unsourced choice must not decide the verdict); sec. 3.10 states that a B3 / B13 ship moves MFT and every
  DPO-derived product by up to (1 - y), so its side-effect record includes every thiol row and the full held-out
  benchmark, and the owner sign-off is required for that reason too.
- v3.2, 2026-11-05 (coordinator's alignment; no other change): omega = 0.419 dex (ln 0.965), the worse of H0's two sets
  (H0 prereg sec. 3.6 / D2, A1 item 3), replaces 0.390 throughout; 0.390 kept only as a named, printed sensitivity
  (likelihood-only: a reweighting / re-solve, no new engine runs). Derived numbers: the re-solve scalings 10^(+-omega)
  2.62x / 0.381x (were 2.45x / 0.407x); a ratio's 2 omega^2 0.351 dex^2, sd 0.593 dex (were 0.304, 0.552); the 0.5x / 2x
  sensitivities 0.210 / 0.838 dex (were 0.195 / 0.780); the adequacy test's |Delta mu| > omega x sqrt(c) thresholds scale
  the same way (0.419 per point, 0.593 per ratio). Unchanged: the tie bands (measurement noise only, never omega) and
  the compute counts (omega's forms need no new runs).
- v3.3, 2026-11-05 (cross-prereg consistency with M1c and M3; no other change): every Martins slope re-expressed on the
  time-averaged pH it is applied on (sec. 4.3.1), with the drift and its source line per source; label-pH forms kept as
  required agreement. P2's s [-0.03, 0.77] -> [-0.02, 0.62] (Martins 2003 Part II / 2003b x 0.79 - 0.81; Martins 2005b x
  1.22 - 1.47, its high-pH arms drift most; Beck unchanged, no drift reported); P3's e N(0.50, 0.19^2) -> N(0.67, 0.18^2)
  (Martins 2005b's k1 x 1.22 - 1.47; the transfer term 0.18 -> 0.145). B3 has no Martins-derived pH quantity. Part (ii)'s
  screen / power check gain the label-band draw sets (within the ~2 h step; the fit path, if it runs, one engine form
  per fitted rival more).
- v3.4, 2026-11-05 (pre-freeze amendment on new data; no other change): Bork 2022 / 2023 (NF alone, pH 5 initial, 130 C,
  second laboratory; dossiers on `lit-2026-11-05-norfuraneol-floors`) read into secs. 0, 2, 4.2, 10, 17. E16 reassessed:
  the pH-flat aerobic law projects >= 13x above Bork at every positive barrier, but hydroxide-only and partial-order
  laws are not separable (the class barrier band alone spans n = 0.37 - 1.07; O2 and Bork's pH drift confound), so **E16
  stays "undecidable on disk"**; Bork is a **band** (for the joint refit's NF sink at acid pH), never a decisive unit.
  Recorded independently: the shipped pH-free NF sink is 3.1 - 4.3x too slow at pH 5 / 130 C. No contrast, unit or
  power-check line is added.
- v3.5, 2026-11-05 (pre-freeze evidence addition for part (i); no other change): Sohn & Ho 1995 (cysteine alone, pH 8,
  NH3 ~14 / ~76 / 82 % at 110 / 150 / 180 C; dossier on `lit-2026-11-05-m1c-bands`) listed as evidence FOR B1 (rank 5,
  measured product yield), with its caveats and the ~30 - 60x deamination-rate conflict with Mizutani 1960; a band source
  (out of every part (i) objective), never decisive. **B1's primary phi_Cys band is unchanged ([0, 1]); Sohn adds
  phi_Cys in [0.76, 1] as a required-agreement form** (a sugar-free pH 8 transfer), budgeted (part (i) ~18 - 64 h serial,
  ~6 - 21 h wall). Register: Sohn vs `AMINE_FATE_BASIS`; Sohn vs Mizutani; Carta 1996's cystine solubility (Hou 7x, not
  10x).
- v3.6, 2026-11-05 (pre-freeze addition; no other change): Bork's supplementary files (branch
  `lit-2026-11-05-bork-si`) leave E16's status unchanged (no headspace, O2, final pH or temperature series; the NF-alone
  series is Fig. 1D replotted, k = 4.91 +- 0.23e-3 /min) and add blocker 4 (NF's self-reaction with its own fragments,
  MGO and glycolaldehyde, roughly 0.012 - 0.16 /min from a 20 mM partner at pH 5: part of Bork's k may be second order and
  not transfer to Mikami's 1.7 mM or pot-level NF; one NF concentration, uncheckable). The joint refit's Bork band is
  scoped "at ~20 mM NF" with a concentration-transfer caveat (sec. 10).
- v3.7, 2026-11-06 (Amendment A1, sec. 20; written after part (i)'s screen and power check, before any number under
  it exists): part (i) closed for shipping (sec. 18.1); one falsification-only re-reading of the stored draws with the
  +-omega corner term replaced by measurement noise (and a lab-repeat offset / shape form for the unbuffered Amadori
  class). Nothing above sec. 18 is edited; part (ii) is untouched by A1.
- v3.8, 2026-11-06 (Amendment A2, sec. 21; written after A1's run, its confirmation fits and an independent
  adversarial review of both (STANDS WITH FIXES), before any number under it exists): the buffered-pot noise term
  derived from Meynier's own statement and read over its defensible range; the kill restated as "B0 + the shipped
  cysteine sinks", the ledger contradiction resolved, one attribution fit; S3's posterior-predictive test; per-signature
  false-kill rates; the commit-message correction; the shipped path written after the numbers. Nothing above sec. 21
  is edited except this line.
- v3.8.1, 2026-11-06 (A2's independent confirmation: CONFIRMED WITH EDITS; sec. 21.9): text and record only. Secs.
  21.2 and 21.6 are annotated in place (marked), not rewritten. Sec. 21.8 is corrected. CF-134 is new. No fit is
  re-run and no shipped output moves.
- v3.8.2, 2026-11-07 (the S3 application, its own PR; sec. 21.10): a record of what shipped, appended; revised the same
  day after the review of #108 with the coordinator's rulings R1 - R3. Nothing above sec. 21.10 is edited except this
  line.

## 20. Amendment A1 (v3.7, 2026-11-06): part (i), falsification only

*Written by the build's author on the coordinator's decision of 2026-11-06, after an independent review of the
noise model. It is committed before any number under it exists, and no amended number is computed until an
independent reviewer who has not seen the screen or the draws confirms this text. One amendment, one re-run; the
result is accepted whatever it is.*

### 20.0 What the author has seen

The author built part (i) and has seen:
- the screen record under v3.6's noise (signature pass fractions, coverage boxes, each rival's best-draw predictions
  and the quantiles printed there);
- the v3.6 power record, and a diagnostic of it (C-Hou's d under v3.6's noise, for three pairs).

**Some amended readings are foreseeable from these predictions.** For example, B0-pinned's Meynier cysteine drifts
(-0.39 to -0.88) stand against a 0.2 bound with a two-reading sd of 0.07. No likelihood, signature, kill or power
number under the noise below has been computed.

### 20.1 Scope

- Part (i) is **closed for shipping** (sec. 18.1): "balance not identified; the measured-trajectory operator stays".
- **Credit outcomes removed.** Outcome rows 1 and 2 of sec. 3.10 ("balance validated" for both classes or for one)
  are removed for part (i) on the data on disk. No rival is favoured, credited, validated or shipped from CF-084's
  controls. The comparison rule of sec. 3.7 items 1 - 2 does not run.
- **The one question left:** does a rival survive CF-084's controls at all? Each rival is either "falsified on
  CF-084's controls" or "not falsified". A rival that is not falsified gains nothing; it stays a candidate for a later
  case on unspent data.
- **Part (ii)** keeps its own text (sec. 4). A1 changes nothing in it. Its pH units must not use `ph_species_discrepancy`
  as built (sec. 20.6).

### 20.2 The error term (replaces sec. 3.3's delta-pH for part (i) only)

- **Primary: per-point measurement noise only.**
  - The sd is the printed replicate SD where one is readable; otherwise max(read error, meter repeatability r =
    0.05 pH). This is sec. 3.3's measurement rule, unchanged.
  - Meynier (C-Mey): the drift is the difference of two readings, so its sd is sqrt(2) x r. It is censored at 0.2, as
    written: P(|d_pred + e| <= 0.2), e ~ N(0, 2 r^2).
  - Buffered pots get nothing more.
  - No model-discrepancy term enters; omega does not appear anywhere in part (i).
- **Required-agreement form: a lab-repeat offset and shape, unbuffered Amadori class only** (C-Hou, C-Hou unforced,
  C-Zhou-ARP, C-Zhou-Cys). Measurement noise as above, plus:
  - a **run-level offset** with sd **s_off = 0.45 pH**, one per run. It is shared by every point of that run: Hou's six
    times share one offset, while each Zhou label is its own run with its own offset.
  - a **per-point shape term** with sd **s_shape = 0.16 pH**, independent per point.
  - So Hou's covariance is diag(meas^2 + s_shape^2) + s_off^2 1 1^T, and each Zhou endpoint has variance meas^2 +
    s_shape^2 + s_off^2.
- **Where the two scales come from (arithmetic only).** They come from one laboratory's repeat of one design: the
  unbuffered AX-ARP 20 mM pot at initial pH 8, 120 C, in Zhou 2024 (Fig. 3b) and Zhou 2025b (Fig. 3d), at 30 / 60 /
  90 / 120 min (sec. 2, CF-060).
  - Readings: 5.13 / 4.32 / 4.19 / 4.08 against 6.04 / 5.05 / 4.66 / 4.50.
  - Differences: 0.91 / 0.73 / 0.47 / 0.42, with mean 0.6325 and sd (n - 1) 0.2296.
  - For a difference of two runs, each run's offset sd is mean / sqrt(2) = 0.447, and the per-point shape sd is
    sd / sqrt(2) = 0.162.
  - The offset scale is a one-pair maximum-likelihood estimate. Its chi-square(1) 90 % interval is roughly 0.23 -
    7 pH, so it is neither an upper bound nor precise. Differences between the papers are part of what it measures,
    which is appropriate: each control is itself a single run in another paper. The shape term includes both papers'
    read noise, so adding meas^2 double-counts about 0.05^2 (negligible). The x0.5 / x2 range does not cover the
    one-pair uncertainty upward.
- **Sensitivities.** s_off and s_shape are scaled together by 0.5 and by 2; these replace omega x 0.5 / x 2. r = 0.02
  is kept, as written. The omega = 0.390 sensitivity is dropped.

### 20.3 The screen's kill rule and signatures (sec. 3.5), noise term only

The signatures, the draws (the stored screen ensembles), the kill rule ("a rival is killed only if no draw reproduces
all four") and the coverage rule are unchanged. Only the noise term inside each signature changes. The readings the
build declared are adopted explicitly, with Z_90 = 1.645 (two-sided 90 %):

| signature | reproduced when (primary) | under the offset / shape form |
|---|---|---|
| sigma-Hou-level | \|pH20_pred - 5.05\| <= Z_90 x meas | Z_90 x sqrt(meas^2 + s_shape^2 + s_off^2) |
| sigma-Hou-shape | the 60 -> 120 min fall within Z_90 x sqrt(2 meas^2) of 0.75 | Z_90 x sqrt(2 meas^2 + 2 s_shape^2); the offset cancels within one run |
| sigma-Mey | every cysteine pot has P(\|d + e\| <= 0.2) > 0.05, e ~ N(0, 2 r^2) | unchanged (a buffered class) |
| sigma-Zhou | the label 6 -> 8 endpoint rise is below B0-pinned's re-printed rise (V6), i.e. toward the measured | unchanged (a sign) |

All four are read on the decisive forced C-Hou run, with the unforced run printed beside it.

**Falsification rule.**
- A rival is **falsified on CF-084's controls** only if the kill holds under the primary term AND under every
  required-agreement form: the offset / shape form at s, s x 0.5 and s x 2; r = 0.02; **r = 0.10 (r x 2; the
  confirmer's edit 3, on sec. 15's rule that an unsourced choice must not decide the verdict: sigma-Mey has no offset
  form)**; and the unforced C-Hou.
- A kill that fails under any form reads "not falsified".
- **Noted openly (from the confirmer).** Under s x 2 the sigma-Hou-level tolerance (Z_90 x sqrt(0.05^2 + 0.32^2 +
  0.90^2) = 1.57 pH) exceeds B0's CF-084 miss (1.30 pH). So in practice a kill can come only from sigma-Mey or
  sigma-Zhou.
- **Confirmation at one fit (sec. 3.8) is kept**, but it is an engine run, so it waits for the coordinator. Until it
  runs, a kill reads "falsified on the screen's draws, confirmation pending".

### 20.4 The power check, re-run as kill power against false-kill rate (replaces sec. 7 for part (i))

- **Truths:** B0-pinned and each rival's best draw. The best draw maximises the summed primary log-likelihood of
  C-Mey (Cys), C-Hou and C-Zhou-ARP under sec. 20.2's measurement noise, on the measured data.
- **Synthetic data:** R = 200 sets per truth, from the truth's stored predictions plus the primary noise.
  - pH points: y = f + e.
  - Meynier: d_obs = d + e, e ~ N(0, 2 r^2), read "within" iff \|d_obs\| <= 0.2.
  - The offset / shape form also gets its own synthetic sets, generated under that form.
- **On each synthetic set** the signatures are re-read against the synthetic values in place of the measured ones.
  sigma-Zhou passes iff the draw's rise minus B0-pinned's has the sign of the synthetic rise minus B0-pinned's. The
  falsification rule of sec. 20.3 is applied to each rival's stored ensemble.
- **Reported per rival k, with Wilson 95 % CIs:**
  - the **false-kill rate**: the fraction of sets from k's own best draw in which k is falsified;
  - the **kill power**: the fraction of sets from each other truth (B0-pinned and every other rival's best draw) in
    which k is falsified; the maximum over truths is printed beside each one.
- **Demotion.** A rival whose maximum kill power does not exceed its false-kill rate is **"unfalsifiable on disk"**.
  Its verdict on the measured data is reported as "not testable here".
- The confirmation fit of sec. 20.3 is not part of the power check.
- **(The confirmer's edit 1.)** Kill power and false-kill rate are reported separately for the primary-noise sets and
  for the offset / shape sets. On each set the full sec. 20.3 rule is applied (all forms). A rival is "unfalsifiable on
  disk" if its maximum kill power does not exceed its false-kill rate (point estimates, as sec. 7) in either family.

### 20.5 No new engine runs; records

- **No new engine run is needed.** Checked on 2026-11-06 from the array shapes only; no value was read. The draws file
  stores, for every draw of every rival, the predicted cooled pH at every scored time of every control pot:
  - C-Mey (Cys) drift: 5 values;
  - C-Hou forced and unforced: 6 each;
  - C-Zhou-ARP and C-Zhou-Cys: 3 each;
  - C-Mey (Gly) for B0, B3 and B13.
  sec. 20.2 - 20.4 need only these and the measured data. The stored delta-pH columns are not read.
- The confirmation fit of sec. 20.3 is the only engine step A1 can call for, and only with the coordinator's go.
- **The re-run writes** `kinetic_core_m1b_part_i_falsification.{json,md}`: the falsification verdict per rival on the
  measured data, and the kill-power table. The v3.6 records (`_screen`, `_power`) are kept unchanged beside it.
- **Run once.** The generator refuses a second run, and refuses to run until this section's confirmation by an
  independent reviewer is recorded in the record's inputs.

### 20.6 `ph_species_discrepancy` (M1.1) before any later pH unit freezes

The function that produced v3.6's delta-pH (`src/kinetic_core/predictive_densities.py`) must be fixed or retired
before M1c's or part (ii)'s pH units freeze. If it is kept:
- co-produced species are scaled together, by their source reaction;
- every scaled amount is capped at the precursor charged;
- the size comes from pH data, never from omega.

### 20.7 Confirmation (2026-11-06)

An independent confirmer returned **CONFIRMED WITH EDITS (text only)**. Its three edits are applied above, in secs.
20.2, 20.3 and 20.4.
- The confirmer disclosed partial un-blinding: it read sec. 18.1 and sec. 20.0's v3.6 figures.
- All three edits make kills harder or clearer, never easier.
- With this, the single falsification run of sec. 20.5 has the coordinator's go, as does the sec. 3.8 confirmation fit
  if a would-be kill appears.

## 21. Amendment A2 (v3.8, 2026-11-06): part (i), the falsification re-examined

*Written by a new session on the coordinator's instruction of 2026-11-06, after an independent adversarial review of
A1's run and its confirmation fits returned **STANDS WITH FIXES** (P ~ 0.80 that B0 cannot reproduce Meynier; P ~ 0.70
on the attribution to the amine ledger). Committed before any number under it exists. Nothing above this section is
edited except the revision-log line (sec. 19). Each step below runs once; its record stands whatever it says.*

### 21.0 What the author has seen

- Everything in sec. 20.0, and A1's records: `kinetic_core_m1b_part_i_falsification` (B0, B2, B3 would-be kills; B1,
  B13 "not testable here"; false-kill rates 0.24 / 0.41 for B0 in the primary / offset-shape families) and the three
  confirmation fits (`..._confirm_{B0,B2,B3}`: kills confirmed). At B0's confirmation fit the Meynier cysteine drifts
  are -0.90 / -1.00 / -0.72 / -0.44 / -0.41 (labels 4.5 - 6.5), the lumped yield is ~2e-8, and the cost is 21.883.
- The review's findings, including one number from its probe: with A1's shape term at s x 2 put on sigma-Mey (sd
  sqrt(2 r^2 + 2 (2 s_shape)^2) = 0.458), **74 of B0's 1280 stored draws pass all four signatures, so B0 would read
  "not falsified"**. The review's probe scripts were read; the output of two of them (the pyrophosphate pKa shift and
  the k_cys_thermal scaling) was not.
- The shipped code of the cysteine sinks (sec. 21.2), read for this amendment; no engine run.
- **Foreseeable:** at B0's confirmation-fit point the smallest drift is -0.41, so any sigma-Mey sd near 0.1 kills that
  point; a sd near 0.46 rescues some screen draws.

### 21.1 The noise term of a buffered pot (replaces sec. 20.2's "Buffered pots get nothing more")

**Stated openly.** Sec. 20.2 gave sigma-Mey only the meter term (sqrt 2 x r = 0.071) after the drifts were known (sec.
20.0), while the unbuffered class got A1's offset / shape term. Put on Meynier at s x 2 (sd 0.458), the shape term
rescues B0 on the screen's draws (sec. 21.0). The buffered-pot term is therefore derived here from measurement-level
arguments, fixed before any number, and the verdict is read over its whole defensible range.

1. **Meynier's own statement bounds it (N-bound).** "No change of more than 0.2 pH unit occurred after heating" covers
   all 20 pots of the paper (glycine, lysine, cysteine and methionine, each at the five labels; meynier1995 dossier).
   - Each reported drift is (true drift + e), e ~ N(0, sigma^2), independent across pots. sigma is everything that is
     not the pot's chemistry: the two readings, the run-to-run spread, the buffer make-up.
   - The most permissive case is a true drift of exactly 0 in every pot. **sigma_max** is the largest sigma with
     P(all 20 |e_i| <= 0.2) >= 0.05, by Monte Carlo (10^6 sets per sigma on a 0.001 grid, seed 20261108), checked
     against the closed form (2 Phi(0.2 / sigma) - 1)^20 = 0.05.
   - A sd above sigma_max is contradicted at the 5 % level by the statement it would be applied to, even if no pot
     drifted at all.
   - **Read-resolution form.** If Meynier read pH to 0.1 unit, a drift up to 0.25 is reported as "0.2". So
     sigma_max,0.25 is computed with edge 0.25, and that form also censors sigma-Mey at 0.25.
2. **Buffer-capacity transfer of A1's shape term (N-beta), per label.** A pH deviation is a proton-count deviation
   divided by the buffer capacity, so the shape term transfers as s_k = m x s_shape x beta_ref / beta_k, m in {0.5, 1,
   2}, and sigma_k = sqrt(2 r^2 + 2 s_k^2).
   - beta_k: Meynier's cysteine pot at label k at t = 0 (0.2 M pyrophosphate, its cysteine and ribose), cooled
     (25 C), at the ionic strength the balance solves, Davies-corrected (sec. 21.6, item 3).
   - beta_ref: Zhou's ARP-alone label-8 pot (the design whose repeat gave s_shape) at B0-pinned's 60-min state, the
     same way. This is the incumbent's state, the only one on disk (declared).
3. **The defensible range D.** sigma-Mey's sd is defensible if it lies in [0.071, sigma_max] at every label (edge 0.2),
   or in [0.071, sigma_max,0.25] with edge 0.25.
   - The verdict is computed at sigma = 0.071, 0.100 and sigma_max (edge 0.2), at sigma_max,0.25 (edge 0.25), and at
     each N-beta set (m = 0.5, 1, 2) that lies in D.
   - N-beta sets outside D and the review's unscaled shape forms (s x 1: sd 0.237; s x 2: sd 0.458) are computed and
     printed, labelled "outside D (contradicted by Meynier's statement)". They do not enter the verdict.
4. **Where it enters.** In every A1 form (sec. 20.3), sigma-Mey's sd is the A2 choice; the A1 form's r and s still
   govern sigma-Hou-level and sigma-Hou-shape as written.
   - A rival is falsified on the screen's draws only if it is killed in every (A1 form x A2 choice in D) cell.
   - **If the kill holds for some choices in D and not for others, the rival is "conflicting / undecided", not
     falsified**, and is carried with B1 and B13 as S3's union (sec. 21.6).
5. The same choices replace sigma-Mey's synthetic noise in sec. 21.4.

### 21.2 Attribution: "B0 + the shipped cysteine sinks"

1. **The ledger contradiction, resolved from chemistry and the cited source (no number needed).**
   - The conflict: `r_cys_thermal`'s docstring (`src/kinetic_core/sulfur.py`, B2.3 charge fix) says the lump's largest
     member, cystine, "RETAINS both of cysteine's titratable centres", and that the pair "had to move TOGETHER". The
     products and `CENTRE_LEDGER` (same commit, 7f65ccaa) carry the carboxyl into CBX, send the nitrogen to FRAG_N and
     declare amine -1 on `AMINE_FATE_BASIS`.
   - The source: Zhai 2023, re-published as Kang 2026 SI Fig. S4 (`KANG_CYS_ANCHOR`), measures the loss of FREE
     cysteine, 10 mM, pH 7, sealed, sugar-free, 100 - 140 C. It names no product.
   - The chemistry: of the lump's named members, cystine keeps both centres of each half (alpha-ammonium pKa ~8.0 /
     8.7, carboxyl ~1.5 / 2.1); thioether self-condensation (lanthionine type) keeps both; a cyclic dipeptide would
     bury both. None destroys the amine while carrying the carboxyl. `AMINE_FATE_BASIS` rests on the amine being
     consumed by Maillard partners into pyrazines, N-heterocycles and melanoidins (Zhou's sugar pots), and the
     measured system has no such partner.
   - **Finding: the docstring is right; the shipped encoding is wrong for this step.** At pH 4.5 - 7 it releases ~1
     proton per cysteine consumed where the step's own chemistry releases ~0. `ch_cys_ox` (autoxidation to cystine,
     "routed exactly as r_cys_thermal") carries the same error. `r_cys_h2s` (H2S + NH3 + pyruvate) is not a ledger
     error: its NH3 is a base that stays titratable unless consumed, and in a sugar pot that is B1's phi_cys question.
2. **The fix, in its own commit.** Flipping the shipped encoding moves the dynamic pH of every cysteine pot, which is a
   shipped output, and this PR moves none (sec. 21.6). So that commit:
   - corrects the docstring and the ledger's basis text to state what the code does and that the chemistry
     contradicts it;
   - adds the chemically right encoding as an opt-in switch: `ph.amine_survival`'s `retained_steps`, whose named steps'
     amine survives as ammonium at fraction 1 in B1's bookkeeping pool (inline coupling only; off by default and
     bit-identical when off), with a test;
   - and the default flip is listed in sec. 21.8 as a never-cased structural choice for S1 / S2.
3. **The claim restated.** A1's kill, where it stands, is a kill of **"B0 + the shipped cysteine sinks"** (their
   encodings and constants), not of B0's balance alone. *[Confirmation edit, sec. 21.9: restated as **"B0 + the shipped
   centre ledger"**. The ledger's carboxyl routing into CBX (CF-134's CO2) is part of what is killed, not only the
   sinks' amine fates.]*
4. **The attribution fit F-A** (one fit, run once).
   - **Rival:** B0, with the shipped encodings.
   - **Free:** the confirmation fit's set (R0's 17 + the lumped yield), with `k_cys_thermal`'s search bound replaced by
     its band, plus mu_h2s, a log10 multiplier on Zheng & Ho's operative `r_cys_h2s` rate (every pH and temperature),
     within its band. The other free constants keep the confirmation fit's bounds; they are not cysteine sinks.
   - **Band of log10 k_cys_thermal** (log10 /min at the engine's 145 C reference; Ea 55.1 kJ/mol as shipped):
     - upper end: Kang's measured total free-cysteine loss at pH 7, k = 8.0e4 exp(-55.1 kJ/mol / RT) /min (the shipped
       anchor), at 145 C. A measurement of the whole of a pool's loss is an upper bound on any one of its channels (S1,
       "bounded");
     - lower end: that x rho, rho = the smallest, over Meynier's five labels, of Zheng & Ho's measured pH law
       k(pH) / k(7) at 140 C (`cysteine_thermolysis_k`). It is the only measured pH law of a thermal cysteine loss (an
       analogue, S1), and the pH-flat constant serves pH 4.5 - 7.
   - **Band of mu_h2s:**
     - upper end: min(0, log10[k_Kang(140 C) / k_Zheng(140 C, pH 7)]). Zheng's constants are flagged as an upper bound
       on free-sulfide release, and their extrapolation cannot exceed Kang's measured total loss at pH 7;
     - lower end: log10 of exp(-(Ea_Zheng,7 - Ea_min) / R x (1/383.15 - 1/413.15)), where Ea_min is the lowest
       measured barrier of the sulfur class (`ea_class_prior`): Zheng's 110 C constant carried to 140 C with that
       barrier instead of its own. This is the 30 C extrapolation past Zheng's data, bounded by the class.
   - **Objective:** B60's 52 rows and prior (the confirmation fit's objective), plus Meynier's five cysteine pots,
     each as a deviance residual sqrt(-2 ln P(|d_k + e| <= 0.2)), e ~ N(0, sigma_max^2) (the most permissive
     defensible choice). *[Confirmation edit, sec. 21.9: inaccurate wording. sigma_max is the largest sd in D. That
     makes it the most permissive pass threshold for a drift beyond the edge, but in the objective it is the weakest
     pull toward Meynier: a smaller sd penalises the same drift more and would push the fit harder toward the five
     pots. F-A is therefore not the fit most favourable to B0 on Meynier. k_cys_thermal sits at its band edge, which
     suggests a harder pull buys little, but no such fit was run.]* The five pots run in-process with BDF at the screen's tolerances (V-solver: BDF and LSODA
     agree to 8e-8 pH).
   - **Start:** B0's converged confirmation-fit point; the two new coordinates at the points of their bands nearest the
     shipped values. M1's optimiser settings and budget, with one continuation of a status-0 segment.
     MAILLARD_WORKERS = 2, memory-gated (>= 35 % free).
   - **Reported:** status; the five drifts; sigma-Mey under each choice of sec. 21.1; the B60 cost (52 rows + prior)
     at F-A and its change from the confirmation fit's 21.883; which band edges are active.
   - **Reading.** "B0 + the sinks meets Meynier" iff every drift has P(|d + e| <= 0.2) > 0.05 at sigma_max. A B60
     cost increase that buys it is reported, not judged here (S2 / S4 would judge it).
   - **Diagnostic only, no verdict:** the drifts at F-A's and at the confirmation fit's points with `retained_steps =
     (r_cys_thermal, ch_cys_ox)`.

### 21.3 S3's test at the converged fit

S3 (STATUS_AND_PLAN sec. 7): a rival is falsified when its converged best fit is inadequate (posterior-predictive p <
0.05 after Holm across the datasets tested). A1's sigma-Mey was a signature screen, not this test.

- **Where.** B0 at F-A (decisive: the rival's best fit with the tested data in its objective) and at its confirmation
  fit (reported). B2 and B3 at their confirmation fits (reported; they have no F-A).
- **Predictive.** The Laplace covariance of the fit's free coordinates (`fitting.laplace` on that fit's objective),
  propagated by a central-difference Jacobian to the control-pot predictions the signatures read (Meynier cysteine 5,
  Hou forced 6, Zhou ARP alone 3).
- **One p per signature.**
  - sigma-Mey: P(all five read within the edge), by Monte Carlo (2 x 10^5 draws of d ~ N(d_hat, J Sigma J^T) plus e):
    the probability of the observed outcome, which is the extreme one (all within).
  - sigma-Hou-level: the two-sided normal p of 5.05 - pH20_hat, sd sqrt(meas^2 + var_pred [+ s_shape^2 + s_off^2]).
  - sigma-Hou-shape: the two-sided normal p of the measured minus predicted 60 -> 120 min fall, sd sqrt(2 meas^2 +
    var_pred [+ 2 s_shape^2]).
  - sigma-Zhou: the two-sided normal p of the measured minus predicted label 6 -> 8 rise, sd sqrt(sd_6^2 + sd_8^2 +
    var_pred [+ 2 (s_shape^2 + s_off^2)]). (The bracketed terms apply under the offset / shape forms.)
- **Holm** across the four (`predictive_densities.holm`, alpha 0.05), in every (A1 form x A2 choice) cell.
- **Reading.** Falsified under S3 iff some signature is rejected after Holm at F-A in every cell with an A2 choice in
  D. If the reading differs across D: "conflicting / undecided".

### 21.4 Per-signature false-kill rates and kill power, every rival

On the stored draws, with A1's synthetic sets (R = 200 per truth and family, seed 20261108):
- per rival and per signature, the fraction of sets from the rival's own best draw in which that signature alone fails
  at every draw (the false-kill rate by signature), and the same from every other truth (kill power by signature),
  with Wilson 95 % CIs;
- the fraction of sets in which Meynier reads all five within;
- under A1's sigma-Mey (sqrt 2 x r) and under sigma_max (the generating noise equals the reading's noise).

### 21.5 The commit-message record

Commit 8eda72ef's message calls the confirmer "the blind confirmer". Sec. 20.7 records a partial un-blinding: it read
sec. 18.1 and sec. 20.0's v3.6 figures. The correct description is **"an independent, partially un-blinded
confirmer"**. History is not rewritten; this line is the correction.

### 21.6 The shipped path

1. Written after the numbers, as sec. 21.8, for a separate PR. **No shipped output moves in this PR.**
   - If B0 + the sinks (restated: B0 + the shipped centre ledger, sec. 21.9) stays falsified: decide between refusing the balance-dependent outputs and a B1 u B13 envelope
     (S3: the union of the carried rivals' envelopes), and state whether B1 and B13 run in the engine today.
   - If undecided: B0, B1 and B13 are carried as S3's union.
   - Either way, list the outputs affected, and list the ledger's default flip (sec. 21.2) as a structural choice for
     S1 / S2.
2. Nothing else in parts (i) or (ii) changes.
3. **Minor (the review's item 6).** `ph_state.buffer_capacity_mmol_per_ph` ignores ionic strength. It gains an
   `ionic_strength_m` argument (default 0: every existing caller unchanged), and A2's betas use the solved ionic
   strength (Davies). CF-084's beta = 5.4 mmol/L per pH at label 4.5 is an I = 0 number; A2's record prints the
   Davies value beside it. Davies beyond ~0.5 M is an approximation (stated).

### 21.7 Records

- `kinetic_core_m1b_part_i_a2.{json,md}`: the noise derivation, the verdicts by choice, the per-signature rates and
  the S3 p-values.
- `kinetic_core_m1b_part_i_a2_attribution_fit.{json,md}`: F-A.
- Each is written once; the generator refuses a second run.

### 21.8 Results and the proposed S3 application

*Appended 2026-11-06 after the three A2 runs; nothing above is edited. Records: `kinetic_core_m1b_part_i_a2`,
`..._a2_attribution_fit`, `..._a2_s3` (the S3 p-values sit in their own record rather than in `_a2` as sec. 21.7
said: a split of the record, not of the method).*

#### 21.8.1 The noise term and the verdict on the stored draws (sec. 21.1)

- **sigma_max = 0.135** (edge 0.2) and **0.168** (edge 0.25); the Monte Carlo agrees with the closed form. So D is
  [0.071, 0.135] (edge 0.2) or [0.071, 0.168] (edge 0.25).
- *[Confirmation edit, sec. 21.9:]* **sigma_max rests on independence across the 20 pots** (sec. 21.1 item 1).
  - A shared error (one meter calibration, one buffer batch, one read-off habit) correlates the e_i. That raises
    P(all within) at a given sigma, so the bound loosens. Fully shared, it would allow sigma up to ~3.
  - The bound is therefore an assumption-laden upper end of D, not a measurement. The S3 verdict does not depend on
    it (sec. 21.8.3): it holds at every sd computed up to 0.458.
  - A bound from the five cysteine pots alone, still assuming independence, is (2 Phi(0.2 / sigma) - 1)^5 = 0.05, so
    **sigma_max,5 = 0.2 / Phi^-1((1 + 0.05^(1/5)) / 2) = 0.265** (checked: 0.05^(1/5) = 0.549, Phi^-1(0.775) = 0.755).
    That still kills (sec. 21.8.3).
- **Davies buffer capacities** of Meynier's cysteine pots at t = 0: 13.2 / 33.4 / 74.5 / 113.6 / 107.8 mmol/L per pH
  at labels 4.5 - 6.5 (I = 0.62 - 1.00 M; at I = 0: 5.6 / 9.4 / 25.1 / 62.6 / 110.1). CF-084's beta = 5.4 at label
  4.5 is the I = 0 reading; Davies gives 13.2 (the review estimated ~12). Davies above ~0.5 M is an approximation.
  beta_ref (Zhou ARP label 8, B0-pinned 60 min) is 2.4. So A1's shape term carried through the buffer ratio is
  0.071 - 0.108 per drift: inside D at every multiplier.
- **Verdict** (draws passing all four, most in any A1 form):

| rival | every choice in D | unscaled s x 1 (0.237, outside D) | unscaled s x 2 (0.458, outside D) | A2 reading |
|---|---|---|---|---|
| B0 | 0 | 0 | 74 | falsified on the screen's draws |
| B2 | 0 | 0 | 15 | falsified on the screen's draws |
| B3 | 0 | 0 | 640 | falsified on the screen's draws |
| B1 | 49 - 95 | 119 | 197 | not falsified (and "not testable here", A1) |
| B13 | 172 - 255 | 273 | 350 | not falsified (and "not testable here", A1) |

- **The verdict does not flip within D.** It flips only at the review's unscaled s x 2 form, which puts sd 0.458 on a
  drift. Under that sd, P(all 20 of Meynier's pots read within 0.2) is ~4e-10 even if no pot drifted at all, so the
  form is contradicted by the statement it is applied to.

#### 21.8.2 Attribution (sec. 21.2)

- **The ledger finding:** the docstring is right; the shipped `r_cys_thermal` and `ch_cys_ox` encodings destroy an
  amine their own chemistry keeps (sec. 21.2 item 1). Fixed as an opt-in switch, commit c26a4a27; the shipped default
  is unchanged. *[Confirmation edit, sec. 21.9:]* the opt-in fix is right for pH <= 7 only. It keeps the amine as
  ammonium (NH4+, pKa 9.25) in B1's pool, whereas cystine's alpha-ammonium pKa are ~8.0 / 8.7, and the ledger ignores
  the thiol centre the oxidation removes (cysteine's thiol, pKa 8.33 at 25 C). Below pH 7 neither matters (both amine
  forms are protonated, the thiol mostly so); at Zhou's label 8 both do.
- **F-A converged** (status 2, 258 evaluations, 0 failed pot runs). Bands: log10 k_cys_thermal [-2.818, -1.980]
  (Kang's total 1.05e-2 /min at 145 C; rho 0.145), mu_h2s [-1.048, -0.083]. Zheng's extrapolated pH-7 rate at 140 C
  (1.05e-2 /min) exceeds Kang's measured total loss (8.6e-3 /min), so the shipped H2S route is above its own bound
  there, as the review suspected.
- **B0 + the sinks cannot meet Meynier:** drifts at F-A -0.50 / -0.64 / -0.55 / -0.37 / -0.36, failing at every choice
  in D. k_cys_thermal sits at its lower band edge and mu_h2s at -0.22. The attempt costs B60 +6.20 (21.883 -> 28.086)
  and degrades Hou's shape (sec. 21.8.3).
- **Diagnostic, no verdict:** with `r_cys_thermal` and `ch_cys_ox` keeping their amine, the drifts are -0.26 / -0.52 /
  -0.51 / -0.36 / -0.35 at F-A and +0.07 / -0.37 / -0.49 / -0.36 / -0.35 at the confirmation fit. The ledger error
  carries the label-4.5 drift.
- *[Confirmation edit, sec. 21.9: the next sentence was a finding; it is a hypothesis.]* **Hypothesis:** most of the
  drift at labels 5 - 6.5 comes from the other cysteine routes (B1's phi_cys question, and the CO2 treatment of CF-134).
  - **What the records show.** A per-step proton budget cannot be computed from the records on disk. The attribution
    record stores, per Meynier pot, the cysteine left and the ammonium the two retained steps made, not the flux
    through each step. Taking the pot's 50 mM start, the records give only the split of the cysteine consumed.
  - At F-A, 10.8 / 20.6 / 32.7 / 37.1 / 38.3 mM of cysteine is consumed at labels 4.5 - 6.5. Of this, 3.2 / 2.4 / 1.4 /
    1.1 / 1.0 mM (29 / 12 / 4 / 3 / 3 %) passes through `r_cys_thermal` + `ch_cys_ox`. At the confirmation fit the
    shares are 71 / 45 / 22 / 16 / 14 %.
  - So at labels 5 - 6.5, 88 - 97 % of the cysteine lost at F-A goes through the other routes. Those routes include the
    three decarboxylations of CF-134, whose CO2 the ledger counts as a full-strength carboxylate. How the drift splits
    between their amine fate and their CO2 is not on disk; it needs a per-step engine run, which this edit does not make.
- **So the kill stands as "B0 + the shipped centre ledger", with the sinks' constants free in their bands.** The
  centre ledger means `CENTRE_LEDGER`'s amine fates and the carboxyl routing into CBX, including CF-134's CO2.
  - P ~ 0.85 (the confirmer's; the author wrote ~0.9) that this structure cannot reproduce Meynier.
  - **P ~ 0.55 on the attribution to the amine fate at release** (the confirmer's; the author wrote 0.85, above the
    review's 0.70). For the attribution: the lumped acid yield fits to 0, and keeping the amine (B1, B13, or the
    ledger's own chemistry) removes much of the drift.
  - Against it: **the glycine pots do drift under B0.** Over the screen's draws the most negative glycine drift per
    draw has median -0.14 and q05 -0.27 (min -0.28; `kinetic_core_m1b_part_i_screen.md`). That pot has no cysteine
    amine, so some of B0's acid comes from elsewhere. B3 and B13 drift by less than 0.006 there.
  - Also against it: the label 5 - 6.5 budget above is not resolved by step, and CF-134's CO2 shares those routes.

#### 21.8.3 S3's test (sec. 21.3)

| fit | cells in D inadequate after Holm | sigma-Mey p (range over D) | other rejections |
|---|---|---|---|
| B0 at F-A (decisive) | 49 / 49 | 5e-5 - 9e-4 | Hou shape (42 / 49; p 4.9e-9 primary) |
| B0 at its confirmation fit | 49 / 49 | 7.6e-4 - 3.0e-3 | Hou shape, Zhou (42 / 49 each) |
| B2 at its confirmation fit | 49 / 49 | 0 - 3e-4 | none |
| B3 at its confirmation fit | 49 / 49 | 0 - 1e-5 | Hou shape (42 / 49) |

- **B0, B2 and B3 are falsified under S3**, on sigma-Mey in every cell, with the Laplace predictive spread included
  (no coordinate unidentified). sigma-Mey stays rejected even under the out-of-D unscaled s x 2 form (p 4e-4 at F-A).
  The reason: S3 asks for all five pots within jointly at the fit, while A1's screen read each pot separately over the
  band-integrated draws.
- *[Confirmation edit, sec. 21.9:]* **The S3 verdict holds over the whole interval [0.071, 0.458]**, not only over D.
  - At F-A, sigma-Mey is rejected after Holm in every A1 form at every sd computed with edge 0.2: 0.071, 0.100, 0.135,
    the N-beta sets, 0.237 and 0.458. Its largest p is 4.7e-4 (0.458, offset / shape form); at edge 0.25 it is 9.1e-4.
    Even Holm's strictest threshold across four tests is 0.0125.
  - **The five-pot bound (0.265) still kills.** It lies between two computed sds: 0.237 (p <= 4.0e-4) and 0.458
    (p <= 4.7e-4). The record keeps the predictive sd per pot, not the 5 x 5 covariance, so p at 0.265 itself cannot be
    recomputed. Under independent pots, P(all within) rises monotonically over [0.071, 0.458] at F-A, which brackets
    p(0.265) at ~3 - 5e-4. That is more than an order of magnitude below 0.0125.

#### 21.8.4 Per-signature false-kill rates (sec. 21.4)

- B0's false-kill rate (0.21 primary / 0.45 offset-shape at either generation sd) comes **entirely from sigma-Zhou**.
  sigma-Mey's false-kill rate is 0.000 for B0 and B2 and <= 0.065 for B3, and its kill power is 0.98 - 1.00.
- Synthetic sets from B0's own best draw never read Meynier "all within" (fraction 0.000). That is because B0's best
  draw drifts -0.4 to -0.9: an unreproducible truth reads "within" at no realistic noise. It is not a defect of the
  check. B1's and B13's own sets read all within 0.14 - 0.34 of the time.
- Hou level and Hou shape kill nothing alone (kill power <= 0.005): the kills are sigma-Mey's and sigma-Zhou's, as sec.
  20.3 foresaw.
- *[Confirmation edit, sec. 21.9:]* **Calibration flag.** These false-kill rates are screen rates under the test's own
  noise model: synthetic sets generated by the same likelihood the test reads. They are not calibrated against
  replicate data, so they bound the test's error only if that noise model is right.
  - One sign that it is not: at B0's own confirmation fit, Hou's shape p is 2e-11 in the primary form (1.5e-13 for B3),
    against 3e-5 - 2e-1 once A1's offset / shape terms are added (1e-2 at s).
  - A p that small at a rival's own best fit suggests the primary form is overconfident. Its rejections are therefore
    read with the offset / shape cells beside them, which is how sec. 21.3 reads them.

#### 21.8.5 The record correction (sec. 21.5)

8eda72ef's "blind confirmer" reads "an independent, partially un-blinded confirmer" (sec. 20.7); recorded in sec. 21.5.

#### 21.8.6 Which shipped path follows (proposed S3 application, for a separate PR; nothing moves here)

- **Status.** The incumbent B0, with the shipped cysteine sinks (restated: **B0 + the shipped centre ledger**, sec.
  21.9), is falsified under S3, and so are B2 and B3. B1 and
  B13 are not falsified. They are S3-carried, since both were named in this case before results with a mechanism of a
  listed evidence kind (B1: Sohn & Ho's measured NH3 yields; B13: Davidek 2006a's acetic yields plus B1). A1 records
  them as "not testable here", so they gain nothing: carried, not credited.
- **S3's text reads two ways here.**
  - The envelope clause: the output is the union of the carried rivals' envelopes, here B1 u B13.
  - The refusal clause: "if the incumbent itself is falsified with no rival credited ... the outputs that depend on it
    are refused". Carried rivals are never credited, so read literally this clause always applies.
  - Both readings agree on one thing: B0's point must not be shown. They differ on whether B1 u B13 may replace it.
    That is for the S3-application review to settle (flagged to the coordinator).
- **Decision: refuse now; the B1 u B13 envelope only through its own reviewed case.**
- **Are B1 and B13 runnable today?**
  - Yes, in the NumPy engine, behind `ph.amine_survival` / `ph.acid_supply`, on the inline coupling only. The shipped
    `fixed_point_nodes` coupling refuses them, and `engine.predict` runs neither by default. They are not in JAX (W3).
  - Their coordinates are bands (phi_cys, phi_arp [0, 1]; y5, y6, y_hi, k_iso), so an envelope is a band-integrated
    ensemble. The screen's 640 draws per rival over ~17 pots took ~6.4 h at up to 4 workers, i.e. roughly 1 - 3 s per
    pot per draw. At 8 workers, a prediction envelope over 2 x 640 draws is ~5 - 10 min per pot.
- **Proposed (the separate PR):**
  1. **Immediately: refuse the balance-dependent outputs** (below) at B0's point. The falsified incumbent's point must
     not be shown, and the envelope is not yet built or reviewed.
  2. **Then: ship the B1 u B13 envelope** for those outputs, built by the separate PR on the inline coupling, its cost
     stated, with held-out coverage before / after and an independent review (S3 / S4's replacement of items 6 and 9).
     Where B1's and B13's envelopes disagree in sign on a lever, the lever's direction is refused (S3).
     *[Confirmation edit, sec. 21.9:]* Only after converged S3 fits of B1 and B13 (A1 recorded them as "not testable
     here", so neither has had S3's test at its own best fit) and after the CO2 ledger audit (CF-134), which B1 and
     B13 inherit.
  3. **The ledger default flip** (`r_cys_thermal`, `ch_cys_ox`: amine kept) is listed for S1 as a never-cased
     structural choice. *[Confirmation edit, sec. 21.9: corrected.]* Under the envelope it is subsumed for these two
     steps: B1 at phi_cys = 1 keeps the amine of the seven steps in `CYSTEINE_AMINE_STEPS`
     (`ph_balance_rivals.py`), and both are among them. It does **not** keep every cysteine amine.
     - Three cysteine-amine steps that `CENTRE_LEDGER` sets to amine -1 lie outside both of B1's sets: `r_ttca_deg`,
       `r_hmf_cys` and `r_cys_actz`. For `r_cys_actz`, the ring nitrogen's destruction is the step's own chemistry.
     - `r_tdp_cys_trap` sits under phi_arp.
     - As a change to the incumbent, the flip needs its own S2 / S4 case; it is not made here.
  4. *[Confirmation edit, sec. 21.9:]* **The CO2 treatment is listed for S1 as a structural choice** (CF-134).
     - The decarboxylations `r_cys_actz`, `r_dpo_ddp` and `r_mgo_cys_strecker` route their CO2 into CBX. So does
       `r_dpo_ddp_h2s` when its switch is on. CBX is titrated as the formic / acetic mixture at full strength (pKa
       ~4.25, `ph_state.TITRATABLE_CENTRES["CBX"]`), where CO2 / HCO3- has pKa 6.35 and degasses.
     - At labels 4.5 - 6 the ledger therefore overstates acid by ~0.5 - 1 H+ per event. B1 and B13 inherit it.
     - Its options: the shipped mixture; carbonic acid at pKa 6.35 retained in solution; carbonic acid with
       degassing (a headspace partition, sealed vs open). None is cased.
- **Outputs affected (first pass; the separate PR computes the exact set).** Only the sulfur lane runs the dynamic
  balance (`engine.py`: the trunk-only, acrylamide and lipid lanes do not). Within it, an output is balance-dependent
  where its sensitivity to the B0 -> B1 / B13 switch exceeds the lane's within-lab term (omega 0.419 dex) for relative
  outputs, or tau for absolute ones. First pass, ordered by expected sensitivity:
  1. the predicted pH trajectory itself (cooled and in situ) in every unbuffered or weakly buffered pot (Amadori /
     ARP pots, Kang-type cysteine pots, pyrophosphate at labels 4.5 - 5.0 where beta is 13 - 33 mmol/L per pH);
  2. the thiols MFT and FFT and the H2S pool in those pots. Every pH-resolved sink and source reads the balance's pH:
     the thiolate loss (Kumazawa), Zheng & Ho's pH-resolved H2S law, TTCA formation and opening, and the Strecker /
     trap steps' pH terms;
  3. ratio outputs across pH ladders (e.g. Meynier's MFT fall from pH 4.5 to 6.5), whose arms drift differently under
     the rivals;
  4. strongly buffered pots (beta >= ~60 mmol/L per pH: 0.2 M pyrophosphate at labels 6.0 - 6.5, phosphate near its
     pKa2). Expected to stay below omega and keep their point predictions; listed so the separate PR checks them, not
     assumed.
- **What remains open:** CF-084 stays open; CF-134 (the CO2 treatment) is open; sec. 20.6 (`ph_species_discrepancy`)
  is unchanged; part (ii) is untouched.

### 21.9 The independent confirmation of A2 and its edits (2026-11-06)

*Appended after an independent confirmer read secs. 21.0 - 21.8 and the records. Verdict: **A2 CONFIRMED WITH
EDITS**. P ~ 0.85 that "B0 + the shipped centre ledger" is falsified under S3; P ~ 0.55 on the attribution to the
amine fate. The edits are applied in place, marked "Confirmation edit, sec. 21.9", and nothing is re-run. No fit and
no shipped output moves.*

1. **Sec. 21.8.2.**
   - "The glycine pots do not drift" was wrong. Under B0 they drift: median -0.14, q05 -0.27.
   - The label 5 - 6.5 attribution is now a hypothesis. The records give the cysteine split, not a per-step proton
     budget.
   - The attribution probability is ~0.55.
2. **CF-134 (new, open; links CF-084).** CO2 from decarboxylations is counted as a full-strength carboxylate.
   - Verified in the code: `sulfur.py` `r_cys_actz`, `r_dpo_ddp` and `r_mgo_cys_strecker` (and `r_dpo_ddp_h2s`)
     route it into CBX. `ph_state.titratable_inventory` titrates CBX at full strength as `ORGANIC_ACID_MIX`.
   - The kill is restated everywhere as "B0 + the shipped centre ledger". The CO2 treatment is an S1 structural choice
     (sec. 21.8.6 item 4). No code changes.
3. **Sec. 21.8.6 item 3.** B1 at phi_cys = 1 does not keep every cysteine amine. `r_ttca_deg`, `r_hmf_cys` and
   `r_cys_actz` lie outside `CYSTEINE_AMINE_STEPS`.
4. **Secs. 21.8.1 / 21.8.3.**
   - sigma_max assumes the pots' errors are independent.
   - The S3 verdict holds over [0.071, 0.458]: at F-A, sigma-Mey's largest p is 4.7e-4, at 0.458.
   - The five-pot bound is 0.265, and it still kills.
5. **Sec. 21.2's objective.** sigma_max is the weakest pull toward Meynier, not "the most permissive" choice for the
   fit.
6. **Sec. 21.8.4.** The false-kill rates are screen rates under the test's own noise model, not calibrated ones.
   Hou's shape p of 2e-11 at B0's own fit suggests the primary form is overconfident.
7. **`fit_m1b_ph_balance.py`.** The S3 step's clip of a stored confirmation optimum to its bounds is now asserted to
   change nothing (`clip_changes_nothing`; unit test). Checked on the three stored records: it changes nothing, so no
   record is regenerated.
8. **The opt-in ledger fix (c26a4a27) is right for pH <= 7 only.** Cystine's amine pKa ~8.0 / 8.7 against NH4+'s 9.25,
   and the thiol centre the ledger ignores, both matter at Zhou's label 8.

**Follow-ups, each its own PR:**
- an S3 application PR that refuses the balance-dependent sulfur-lane outputs;
- a sec. 7 clause on how a refusal is lifted;
- the B1 u B13 envelope, only after converged S3 fits of B1 and B13 and the CO2 ledger audit (CF-134);
- a case for flipping the ledger default.

### 21.10 The S3 application (2026-11-07; its own PR, `s3-ph-balance-refusal`, #108)

*Appended by the session that built the application, revised after the independent review of #108 (MERGE WITH FIXES)
and the coordinator's rulings R1 - R3 of the same day. Nothing above is edited. The sec. 7 lift clause it implements is
"pending independent confirmation".*

**R1. Which term judges which output (coordinator's ruling).** Each output's spread over B0 and the carried rivals is
compared with the uncertainty the output is SHIPPED with: tau = 0.964 dex (H0's tau_ABS(all); sec. 21.8.6, confirmed)
for an absolute concentration shipped as a level, and omega = 0.419 dex for a within-laboratory relative output (a
ratio between pots, a within-study contrast, a direction), with S3's sign rule. The argument: S3 refuses an output that
DEPENDS on the falsified structure, and its tolerance sentence already assigns the within-laboratory term to relative
outputs and the between-laboratory term to absolute ones. A level is shown with an interval that carries tau; a rival
spread inside it moves the level by less than the uncertainty it is shown with, so the output does not depend on the
falsified choice beyond what it already shows. A ratio is shown with omega in its interval (H0's ratio interval), so
omega judges it. The first build used omega for every output (the refusal clause's parenthesis read alone); R1 reads
the parenthesis with the tolerance sentence. omega on absolute levels is reported as a sensitivity.

**R2. pH levels.** The cooled and in-situ pH the balance solves are levels, judged in every pot by the pH-level term,
A1's run offset and per-point shape in quadrature, sqrt(0.45^2 + 0.16^2) = 0.478 pH (sec. 20.2). The 0.135 pH of sec.
21.8.1 is the sd of a DRIFT (a difference of two readings) and rests on the independence of Meynier's twenty pots (sec.
21.9: "assumption-laden"), so it never judges a level; the five-pot bound sigma_max,5 = 0.265 pH is reported as a
sensitivity.

**R3. The rivals' representation.** Each carried rival is represented by its screen draws that pass all four
signatures in at least one (A1 form x A2 choice in D) cell of sec. 21.1 (B1: 95 of 640; B13: 255 of 640,
the counts sec. 21.8.1 reports), through the passing draw that holds the minimum and the one that holds the maximum of
each free coordinate, plus the rival's best draw (B1: 8 points over 4 coordinates; B13: 15 over 7;
`results/validation/ph_balance_s3_rival_points.json`). An output is refused on the largest spread over B0 and all of
them. Why the extremes represent the band: the rival coordinates enter the charge balance as supplies that move the pH
one way each (more surviving ammonium, phi_cys / phi_arp, raises it; more lumped or acetic acid, acid_yield / y / k_iso,
lowers it; a higher ARP pKa holds more base), and every pH-resolved rate reads that pH, so an output is expected to be
monotone in each coordinate across the passing box and its extremes sit at the box's edges. **Checked** on the record's
pots (8 pots, 24 further passing draws per rival, seed 20261107): **the premise fails.** 23 of 58 outputs' spreads widen beyond B0 + the extremes when further passing draws are added (largest widening 2.05 pH units, the cooled pH of cysteine + ribose + glycine in 0.2 M acetate at pH 7.5, 0.71 -> 2.76), though no verdict changed in that sample: the coordinates act jointly, so a corner of the box can lie beyond every one-coordinate extreme. Run against EVERY passing draw (B1 95, B13 255) on the same pots (`draws_needed` in the record): the extremes reach within 0.05 of the all-draws spread in 0.76 of the 58 outputs and give the all-draws verdict in 57 of 58 (the one miss: the MFT disulfide of the ARP + cysteine pot in 0.2 M phosphate, 0.94 against 0.998 dex at tau 0.964, a level H0 refuses anyway). Random subsets need more than 80 draws per rival to reach 0.95 (0.34 / 0.43 / 0.54 / 0.68 / 0.83 at 5 / 10 / 20 / 40 / 80; verdicts agree in 0.92 - 0.98). **So the band is bounded only by all 350 passing draws per pot** (~0.4 s each: a median of ~140 s per pot in the NumPy engine, against ~10 s for the 23 extremes). The engine ships the extremes as R3 rules; that they under-cover the band, and the cost of using every passing draw, are open for the coordinator (the envelope case of sec. 21.8.6 item 2 would carry the draws anyway).

**Fix 1 (the review's main issue).** A level refused under R1 stays in the arm, withheld by name as H0 withholds a thiol
level (`absolute_refused_targets`, with `balance_refused_targets` for the reason; no calibration rescues it), so every
ratio, contrast and direction is judged by the ratio rule (omega and the sign) and never by its arms' levels. compare's
state vectors drop the refused species (fix 2).

**Fix 4.** A rival integration is stopped after a deterministic budget of right-hand-side and Jacobian evaluations
(`ph_balance_rivals.RIVAL_RHS_BUDGET` = 200 000; every finished rival point needs at most ~38 000 on the slowest pots
probed), not a wall clock; a pot whose rival exceeds it is refused with that reason. In the record: five arms are refused for want of a rival: Chen 2000 / 2002's four pH 5 ribose + cysteine arms (180 C, 2 h) and Hofmann's dry-heated 180 C / 6 min ribose + cysteine (HOF-03).

**What the record shows** (`results/validation/ph_balance_s3_sensitivity.md`, 809 recipes: 26 sulfur-lane panel pots,
165 arms of 65 directional claims, 618 grid recipes incl. 88 clamped controls and 2 lipid co-integrations):

| class (balance-solved pots) | refused / decided | at the sensitivity term |
|---|---|---|
| absolute concentrations (tau 0.964 dex) | 1155 / 3522 (33 %) | omega 0.419: 1885 (54 %) |
| cooled + in-situ pH (0.478 pH) | 671 / 1396 (48 %) | 0.265 pH: 813 (58 %) |
| ratios, pH lever (omega, sign) | 1201 / 2603 (382 sign flips) | - |
| ratios, buffer lever | 1697 / 2499 (585) | - |
| ratios, dose / sugar / thermal | 285 / 374, 129 / 607, 337 / 1759 | - |
| within-study panel pairs | 57 / 102 (19) | - |
| directional claims (sulfur lane) | 27 / 45 (12 direction flips) | - |
| clamped controls | 0 / 882 | - |

Under fix 1 SUG-04 comes back (its step stands under the ratio rule); Z23C-GSH-01 stays refused (its step's spread is 1.07 dex, beyond omega).

**Not refused.** A clamped (held, measured) pH: the balance is not solved (the measured-trajectory path); the trunk,
acrylamide and lipid lanes. The WUR blind-test predictions of 2026-10-20 ran with `BufferSpec(kind="clamped")` (pH 6.4
held, their assumption A4): they did not depend on the balance and are not re-scored.

**Cost.** One recipe at all 25 points timed at 11.9 s projected 0.67 h for the 809 recipes at 4 workers; the run took 2.03 h (directional arms with long programmes run 70 - 900 s of rivals each). The live rule costs a median 9.6 s of rival integrations per balance-solved pot (23 points, ~0.4 s each; max 910 s, a Wang 2012 sweep arm), the directional generator about an hour at 4 workers.
