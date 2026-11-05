# M3 pre-registration: the Strecker-aldehyde layer. One shared ARP-Phe decay and two routes (phase 1), then the layer fed from the network (phase 2)

**Status: v2 DRAFT; freezes after its power check, before the keyed Figs. 2 - 3 are committed.** Nothing has been
built, screened, fitted or scored. No engine and no layer was run to write this file.

- **What v2 is.** It answers an independent design review of v1 (a66fe933), which found MORE FIXES and asked for no
  redesign. The coordinator's decisions are applied (sec. 21).
- **What must happen before phase 1 freezes,** in this order:
  - the build and its V-checks;
  - the U3 synthetic power check on the expected design (sec. 12.1);
  - the mechanism screen (sec. 8) on every unit except U3 and sigma-trap;
  - the full power check (sec. 12);
  - the answers to any re-review.
- **The keyed Figs. 2 - 3 come after the freeze.** Two blind readers have already digitised them. Their files sit in the
  coordinator's scratchpad. The author has not opened them, nor the PDF's Figs. 2 - 3. Those rows are committed only on
  top of the frozen v2 commit (sec. 7).
- **After the screen, the only change allowed is demoting a contrast** to "undecidable on disk". No band, unit, form or
  rule may change after it (sec. 13).
- **Phase 2 is a pre-declared conditional** (sec. 14). Its units and bands are frozen by an amendment written after
  phase 1's verdict and before any phase 2 run.

*Plan step M3 ([STATUS_AND_PLAN](../../docs/guides/STATUS_AND_PLAN.md) sec. 5, W2, row M3; sec. 7 items 1 - 9; the
re-sequenced plan on branch `plan-reseq` has the same row). Case:
[M3_strecker_aldehydes.md](../../docs/reference/cases/M3_strecker_aldehydes.md) (items 1 - 5, 2026-11-04; not edited
here). Earlier records: EX-042, EX-049, EX-051 ([experiment_register.yml](../../docs/reference/experiment_register.yml)),
[kinetic_core_strecker_forced_input_prereg.md](kinetic_core_strecker_forced_input_prereg.md) and its review, the
feasibility and blind-Fig. 1 records, and B22 / B22b / B24 / B29 (their ship rules here). The structural template is
the M1b prereg v3.1 (branch worktree `agent-a1611d73ae9435878`, unmerged), and this design follows its rules 1 - 16 as
the brief numbers them.*

- v1 2026-11-05 (a66fe933); v2 the same day. Branch `m3-strecker-prereg` from origin/main at cc5cbb55.
  `STATUS_AND_PLAN.md`, the registers, the case files and `data/` are not edited. What this design needs from them is
  listed in sec. 19.
- Every computation runs through `./scripts/docker_maillard.sh run`, `nice`d, on <= 8 cores.
- "(arithmetic)" marks a number computed here from printed or keyed values. **"(my arithmetic)"** marks v2's
  acid-load and O2-budget numbers. Their script, `m3_acid_load.py`, is in the session scratchpad and was run through the
  wrapper. No such number sets a band or a start (sec. 5).

## 0. Seen by the author, and which pots have been seen under the engine or the layer (rule 16)

**The expectations in sec. 16 are not blind.** Read while writing:

- **Layer outputs (EX-049, EX-051), read in full.** These cover the glucose + Phe pot shapes on hofmann2000 Fig. 2
  (DC_ARP optima: 10-min point -0.26 dex, R30 0.38 - 0.49, ARP share of PA at 120 min 45 - 65 % / 0.84 at the search
  optimum) and hofmann2000b Fig. 1, both curves. G5' passes blind and G7 fails blind: +0.5 dex at 10 min and -0.4 to
  -0.5 dex at 500 min, with k_d 0.058 - 0.069 /min.
- **Engine outputs (whole-model refusals), from the ship-rule records.** B29 on hofmann2000b Table 2, both pots (model
  air / argon 2.19 vs 9.2 for the fed ARP, 2.76 vs 3.5 for glucose). B22 on Pan 2025 and Deng 2022 (methional at 120 min
  +1.93 dex). B22b on Deng 2022's shape (the peak at 30 min against 120). B24 on hofmann2000 Tables 7 / 9 (proline).
- **Measured values, no model number.**
  - hofmann2000b Tables 2 - 3 and Fig. 1 (`strecker_hofmann2000b_fig1.yml`): fed ARP PA 0.99 / 2.61 / 4.21 / 6.21 /
    8.49 / 9.20 / 10.00 and PAA 0.30 / 0.80 / 1.81 / 3.21 / 6.01 / 7.20 / 8.99 umol/mmol at 10 - 500 min;
  - hofmann2000 Tables 3 - 4 and the Fig. 2 values quoted in the case;
  - Deng 2021 Fig. 3 / 6 and Table 1, and Deng 2022 Table 1 (dossiers);
  - Martins 2003 Part II Table 3 and Martins 2003b Figs. 6 - 7 (acids, pH drop);
  - Davidek 2006a's acetic-acid yields (dossier);
  - Weenen 2001 Tables II - III.
- **hofmann2000b Figs. 2 - 3: only the dossier's text has been seen**, never the figures. That text gives the dominance
  ordering (1-DG main under argon; glucosone main under air + Cu; 1-DG and 3-DG "only low yields" under Cu; glucosone
  "very low" under argon) and an end point ("Phe liberated ~55 % after 320 min" in both atmospheres). **So the U3-Phe
  end point and the U3 dominance ordering are SEEN; only the time courses and the yields are blind** (M10). The readers'
  files and the rasters were not opened. No model has been computed on these pots.

**Pots seen under the engine or the layer** (no "blind" claim is made for any of them):

| pot | seen under | record |
|---|---|---|
| hofmann2000 Fig. 2 (glucose + Phe, PA / PAA 10 - 500 min) | the forced layer (EX-049 primary and early-step) | feasibility records |
| hofmann2000b Fig. 1, squares and circles | the forced layer (EX-051 blind test) | blind-Fig. 1 record |
| hofmann2000b Table 2, both pots, three atmospheres | engine (B29 fit rows: the PA ratios) and the layer (fed ARP rows were exact calibrations in EX-049) | B29 ship rule; feasibility |
| hofmann2000 Tables 3 - 4 (fed donors) | the layer (exact calibrations) | feasibility |
| Deng 2022 (methional, ARP vs sugar, 120 C) | engine (B22 T3, B22b T1) | B22 / B22b ship rules |
| Pan 2025 (methional rates) | engine (B22 T1) | B22 ship rule |
| hofmann2000 Tables 7, 9 (proline) | engine (B24) | B24 ship rule (out of scope here) |
| hofmann2000b Figs. 2 - 3 (trapped osones) | **none**; the end point and the ordering are seen as text | -- |
| Weenen 2001 Table II, Martins 2003 / 2003b / 2005b, Deng 2021, Davidek 2006a | none (measured values read, no model number) | dossiers |

## 1. Question, claim and phases

### 1.1 The claim

Phenylacetaldehyde (PA) and phenylacetic acid (PAA) form by two routes:

- **Route A (dicarbonyl).** Donor D + free Phe -> PA / PAA, for D in {MGO, GO, 3-DG, glucosone, 1-DG, the fragment
  osones}. This is the EX-049 rate law, with the omega acid channel by class (arm A').
- **Route B (Amadori, oxidative).** ARP-Phe + O2 -> PA / PAA with no free dicarbonyl and no free Phe (hofmann2000b
  Fig. 4).

ARP-Phe decays first order with ONE total constant k_d (100 C, label pH 7.0, 0.5 M phosphate, air), split into exits
whose shares are properties of the compound, not of the pot. **The same k_d and shares hold in the fed ARP-Phe pot and
in the glucose + Phe pot.** So route B's rate in the glucose pot follows from the fed pot with no free parameter.

### 1.2 Phases

- **Phase 1.** This is a standalone forced-input layer, not the engine, extending EX-049's module. The glucose pot is
  forced by its measured courses. The fed ARP pots (closed vial, reflux, in-situ trapped) are solved from the ARP's
  exits. Phase 1 is **single-laboratory** (DFA Garching), and every verdict sentence says so. **A phase-1 layer ships for
  relative within-laboratory predictions only**; absolute Strecker levels for a new laboratory are refused
  (coordinator's decision on v1's owner item 2).
- **Phase 2 (sec. 14, conditional).** It runs only after outcome row **1, 1b, 2 or 3** (sec. 10.3). The phase-1 structure
  is fed by the network. It shares the pyrazine steps' donors and aminoketones, so it moves shipped pyrazines, and
  **the owner signs off at ship time.**

## 2. The evidence: additions and corrections to the case

| # | source | what it adds | consequence |
|---|---|---|---|
| 1 | Martins 2003 Part II Table 3 (fed DFG 10 mM, 0.1 M phosphate, 100 C, Model 2): **k1 (DFG -> E1, the 1,2-enaminol that leads to 3-DG through k4) 0.0057, k2 (DFG -> E2, the 2,3-enaminol that leads to 1-DG) 0.0156, k3 (DFG -> C_n, C5 / C6 fragments; MG directly in Model 1) 0.0155 /min at pH 6.8** (sum 0.0368); sum 0.0047 at "pH 5.5" | the exit shares: 3-DG route 0.15, 1-DG route 0.42, fragments 0.42. The pH slope is 0.69 dec / pH on initial label, but **the "5.5" arm drifted 0.8 - 0.9 units** (Martins 2003b, text) while the 6.8 arm drifted 0.19 | a **band source** for k_d and the shares' start (sec. 4.2). The slope re-derived on the mean label-scale pH of each arm is **0.54 - 0.56 dec / pH** (arithmetic). R-D2X is reported only |
| 2 | Deng 2021 (Met-ARP 0.2 M, 100 C): 0.0117 /min unbuffered at pH 7.5 vs 0.0048 at pH 5.5 (**0.19 - 0.20 dec / pH at 100 C**); 0.020 - 0.022 in 0.1 M phosphate at pH 7.5; Ea "75 +- 10 kJ/mol"; 88 % Met returned | k_d's buffer and pH transfer, and the exits' Ea | a **band source**. **The 1.8x phosphate factor carries a pH confound**: the unbuffered pot's pH is not reported and falls fastest early, so part of the 1.8x may be acidification, not catalysis. R-Deng is reported only |
| 3 | **Correction to the plan row's D2 sentence (downgraded in v2).** At <= 0.1 M phosphate, no measured amino-acid ARP decay at 100 C reaches 0.058 - 0.069 /min. At Hofmann's 0.5 M, if phosphate catalysis is linear (EX-058's form, measured on a different step; declared transfer here), Deng's two points give ~0.057 /min and Martins' 0.037 x 2.8 gives ~0.10 | **EX-049's in-pot values are consistent with** the measured decays at 0.5 M under a linear phosphate law. That law rests on Deng's 1.8x, which item 2's caveat weakens | D2 cannot by itself reject R3 |
| 4 | **hofmann2000b Table 3, read tightly.** The glucose pot's ARP-Phe falls 48 -> 32 umol over 30 -> 60 min. Because formation >= 0, ARP(60) >= ARP(30) e^(-30 k), so k_d,glc >= ln(48/32)/30 = 0.0135 /min; **0.0126 - 0.0144 with integer rounding** (arithmetic). The dossier's 0.0062 (30 -> 300 min) is weaker | a lower bound on the in-pot decay that needs no calibration (the scale cal_ARP cancels). **It assumes one vessel sampled over time and analyser noise within rounding**; sec. 6 integrates a 2 - 10 % analyser CV and a between-vessel form, and both weaken it | sigma-ARP and U-ARPglc |
| 5 | EX-049 put the fed ARP's argon PA on an **anoxic direct route** (k_a0, psi0) | that layer had no route A in the fed pot, so it could not express R2 | this layer solves the fed pot from the ARP's exits |
| 6 | EX-049 calibrated psi on Table 2's air excess | the fed ARP's air PA / PAA was in-sample | psi is integrated over [0, 1] |
| 7 | Granvogl 2006 (DFA, same [13C2]-PA SIDA): PA within-run RSD 1.7 - 5.9 % | replicate noise | s_rep = 0.025 dex |
| 8 | Granvogl 2006 model B (MGO + Phe 5 mM, 170 C, 30 min): PA 72.8 mmol/mol vs Hofmann's fed MGO 11.1 umol/mmol | the same laboratory's spread of efficiency across conditions, x6.6 (0.82 dex) | the upper end of s_DG1. **It is generous to R2**: 170 C, 5 mM and a short time are far from Hofmann's 100 C pots |
| 9 | Weenen 2001 Table II, pH 7.0: ARP / (Glc + Phe) PA x4.2 (Hofmann Table 2 air x3.9) | a second laboratory's scale-free ratio | phase 2 only (no forcing for Weenen's sugar pot) |
| 10 | Deng 2022 | no ARP and no keyed donor in its sugar pot | phase 2 only |
| 11 | **Acids from ARP decay.** Martins 2003b Fig. 6: acetic + formic acid ~6.8 mM from 10 mM DFG (100 C, pH 6.8, 180 min); its measured pH drop of 0.19 in 0.1 M phosphate implies **~1.08 mol acid per mol DFG (my arithmetic)**, as Martins states ("other acids might be formed"). Davidek 2006a: 1-DG gives 0.56 - 0.84 mol acetic acid per mol (pH 6 - 8, 120 C, no amine); DFG gives 0.28 / 0.56 mol acetic at pH 6 / 8 (Davidek 2005 via 2006a) | the acid yield per mol ARP decayed: **y_acid in [0.4, 1.1]** | the in-situ pH estimate (sec. 3.4) |
| 12 | Voigt 2010 (fed 1-DG + Ala, 37 / 50 C, pH 7.4): acetaldehyde x1.5 - 1.6 aerated (the deaerated arm carries DTPA) | O2 can raise a 1-DG Strecker yield | the R2 form "s_DG1 x [1.5, 1.6] in air" (sec. 4.3) |

## 3. The layer (phase 1): structure common to every rival

The module follows `docs/guides/A2_MODULE_CONTRACT.md` sec. 3, in a new generator,
`scripts/generators/m3_strecker_layer.py`. The frozen EX-049 / EX-051 generator is not edited. Nothing in `src/` or
`data/` changes.

### 3.1 Route A, carried from EX-049 (reviewed in EX-042; not re-derived)

- **Rate laws.** The Strecker flux is J_D = k_D [D][Phe]. PA gets (1 - phi_D) J_D, minus k_loss [PA]. PAA gets
  (phi_D + omega_class(D) O2) J_D. The donor sink is -kappa_D [D][Phe], and Phe loses nu kappa_D [D][Phe].
- **Exact inversion per draw.** k_D, phi_D, omega_C23 and omega_C6 are inverted exactly from hofmann2000 Tables 3 - 4.
- **Glucose pot.** It is forced by Fig. 3's donors, by hofmann2000b Table 3's Phe and by ARP-Phe, each course with its
  calibration coordinate.
- **CF-070.** The Table 3 Phe rows of hofmann2000 stay out. Constant Phe is a required-agreement form.

### 3.2 The ARP-Phe exits (shared by the fed pots and the glucose pot)

At 100 C and the pot's in-situ pH (sec. 3.4), with O2 the pot's oxygen state (sec. 3.3):

| exit | product | rate | Phe released |
|---|---|---|---|
| e1 | 1-DG (DG1) | a1 k_anox | yes |
| e3 | 3-DG | a3 k_anox | yes |
| eF | fragments, then MGO (fraction m, yield 1) or inert C | aF k_anox | fraction r_F |
| eG | glucosone | g_G k_ox O2 c_Cu | yes |
| eFox | **oxidative fragments**: GO, MGO or a C4 / C5 osone (lumped as "FOX", a donor whose Strecker share s_FOX is drawn within the **measured** DFA donor range, glucosone to MGO) | g_F k_ox O2 c_Cu | yes |
| eB | **route B**: PA (1 - psi) or PAA (psi) | g_B k_ox O2 c_Cu | no |

- **The constants.** k_d (air) = k_anox + k_ox. The O2 share is sigma_ox = k_ox / k_d in [0, 1). The anoxic shares
  (a1, a3, aF) and the O2 split (g_G, g_F, g_B) each sum to 1.
- **The fed pot.** ARP(t) follows dARP/dt = -k_d,eff(t) ARP, where k_d,eff sums the exits at the pot's O2 and pH.
  Released Phe, plus any free Phe in the ~10 % impurity, feeds route A on every ARP-made donor. That is a small ODE
  per pot.
- **The glucose pot.** ARP is forced. Exit fluxes are (exit rate at k_d,glc) x ARP_forced(t). 1-DG and FOX are
  computed from ARP. The measured 3-DG, glucosone and MGO already contain their ARP-derived parts (no double count).
- **The trapped pot (hofmann2000b Figs. 2 - 3).** ARP-Phe 0.5 M at 85 C, with 1,2-diaminobenzene 0.6 M, under argon or
  air + Cu 25 mM.
  - Exits at 85 C: k_e(85) = k_e(100) exp(-Ea/R (1/358.15 - 1/373.15)), 0.363 at 75 kJ/mol.
  - **Each species is trapped with an efficiency set by competition** (M2). eps_X = k_tr / (k_tr + k_loss,X(85)),
    where k_tr, the pseudo-first-order trap rate at 0.6 M, is in [0.1, 10] /min (unsourced, numerical ends: a verdict
    that depends on them does not ship). For 1-DG, k_loss,X is k0_DG1 + kappa_DG1 [Phe]; for 3-DG and glucosone it is
    the route-A sink.
  - The trap gives cumulative trapped osone. Route A is suppressed by what is trapped. Phe liberated is the sum of the
    releasing exits.

### 3.3 Oxygen: states and a finite budget (M1)

- **Closed air vial:** O2 = 1 in the primary form, an unlimited supply.
- **Finite-budget forms (required agreement).** Headspace air is 5, 20 or 100 mL, at 8.6 umol O2 per mL (25 C, my
  arithmetic), plus ~2.6 umol dissolved in 10 mL. Each oxidative event (eG, eFox, eB and each omega acid event)
  consumes nu_O2 in {0.5, 1} mol O2. O2 falls as it is consumed, and the oxidative rates scale with O2 / O2(0). The
  100 mL form covers the "flushing with oxygen" reading of Table 2 (pure O2 in ~21 mL holds the same as 100 mL air).
- **Argon vial:** the residual o_Ar is in {0, 0.1} (form).
- **Reflux:** O2_reflux in [0, 1], one value shared by every reflux pot. A per-pot O2 is a form for R3 verdicts.
- **Cu:** c_Cu at 5 mM is set by exact inversion on Table 2's air + Cu fed-ARP total, **bounded in [1, 25]** (Cu never
  slows oxidation; the upper end is numerical). A draw outside the bound is infeasible. c_Cu at 25 mM is either
  saturated (= c_Cu5) or linear (1 + 5 (c_Cu5 - 1)).

### 3.4 In-situ pH (B1): an acid load against the buffer, per pot

- **The buffer model.** Phosphate is described by Henderson-Hasselbalch at pKa2' = 6.79 (HPO4^2- fraction 0.62 at
  pH 7.0 in 0.5 M; my arithmetic). The pH falls as the acid formed titrates HPO4^2-.
- **The acid.** Acid formed = y_acid x (ARP decayed), with y_acid in [0.4, 1.1] (sec. 2 item 11). In the trapped pot,
  only osone that is not trapped yields acid: (1 - eps_X) x y_acid per osone, plus the fragment exits.
- **The glucose pot.** Acid = y_acid x (Phe consumed, forced; the sugar-degradation proxy). Under the CF-070 form
  (constant Phe), acid comes from ARP turnover alone.
- **The scale (rule 13).** Every Delta pH is computed on the 25 C buffer equation, from the pot's own label pH 7.0. The
  label-to-in-situ shift is taken to cancel in the Delta (same buffer, same temperature within a pot; declared). The
  band sources' slopes are also on the label scale.
- **Its effect on k_d.** k_d(t) = k_d,7 x 10^(s_pH Delta pH(t)), with s_pH in [0.20, 0.54] dec / pH: Deng's 100 C slope
  and Martins' drift-corrected slope.
  - This is **R4's structure**.
  - For R0-pinned, R1, R2 and R3 it is the required-agreement form **"pH drift on"**, applied in the trapped pot. In the
    0.1 M pots, R1's drift effect is printed beside each unit as a diagnostic.

**The acid-load numbers per pot (my arithmetic; y_acid 0.4 - 1.1):**

| pot | phosphate | acid load at the scored times | pH drop (label scale) | lower k_d late, dex (s_pH 0.20 - 0.54; 0.69 on initial label) |
|---|---|---|---|---|
| fed ARP 0.1 M, closed vial / reflux | 0.5 M | 120 min (~81 % decayed at 0.014 /min): 32 - 89 mM; 500 min (~100 %): 40 - 110 mM | 0.12 - 0.31 at 120 min; 0.14 - 0.39 at 500 min | 0.03 - 0.21 (0.27) |
| trapped ARP 0.5 M (Figs. 2 - 3), if no acid were suppressed by trapping | 0.5 M (as the footnotes print) | 320 min (50 - 55 % decayed): 100 - 300 mM | 0.35 (y 0.4) to 1.1 - 1.3 (y 1.0); at y 1.1 and 55 % HPO4^2- is nearly exhausted (pH <= ~5) | 0.07 - 0.70 (0.90); trapping lowers this |
| trapped pot under the Methods' "0.5 mmol/L" phosphate reading | 0.0005 M | -- | effectively unbuffered (exhausted at 0.06 % acid per ARP) | not computable: a form "unbuffered trapped pot" is declared, and U3 cannot credit under it |
| glucose + Phe 0.1 / 0.1 M | 0.5 M | 120 min: 20 - 54 mM; 300 min: 24 - 66 mM | 0.07 - 0.19; 0.09 - 0.23 | 0.01 - 0.12 (0.16) |
| check: Martins fed DFG 10 mM | 0.1 M | 6.8 mM measured acids | predicted 0.12 vs measured 0.19 | (the gap gives the 1.08 mol/mol implied yield) |

The reviewer's numbers (trapped pot pH ~5.9 - 6.3; 0.1 M ARP pots 0.1 - 0.2) lie inside these ranges. Mine are wider at
both ends, because y_acid spans 0.4 - 1.1.

**B1's consequence: the trapped pot reads k_d at a falling pH.** Its decay can differ from the 100 C pots' by
acidification as well as by temperature and order.
- **The fed block has two conditioning forms** (sec. 6). The full form uses U3. The "100 C fed only" form excludes U3.
- **Required agreement holds between them.** If they disagree, the transfer contrast reads **"confounded by
  acidification"** (row 5).

## 4. Rivals, coordinates, bands, starts

### 4.1 Coordinates every rival carries (integrated over their bands; never fitted to a scored unit)

| coordinate | band | source |
|---|---|---|
| log10 c_MGO, nu, log10 k_loss, off_kappa_GO / G / TDG (with floors) | EX-049's | forced prereg sec. 2; EX-042 review. The +0.4 kappa ceiling is CF-070-dependent; +0.5 is a form |
| cal_MGO, cal_GO, cal_TDG, cal_G, cal_ARP | +-0.30 / 0.40 / 0.20 / 0.20 / 0.20 dex | EX-049 |
| O2_reflux; o_Ar; headspace and nu_O2 | [0, 1] shared; forms {0, 0.1}; forms {inf, 5, 20, 100 mL} x {0.5, 1} | sec. 3.3 |
| arp0; imp_Phe | [90, 100] mM; [0, 1] | purity "about 90 %"; definition |
| anoxic shares (a1, a3, aF); O2 split (g_G, g_F, g_B) | Dirichlet(1, 1, 1) each | definition. Start 1 uses Martins' shares |
| r_F; m | [0, 1]; [0, 1] | definition |
| log10 k_d,7 (at label pH 7.0) | **[-1.95, -0.88] /min** (0.0113 - 0.132); floor form log10 0.005 = -2.30 | sec. 4.2 |
| s_pH | [0.20, 0.54] dec / pH | Deng 100 C; Martins drift-corrected (sec. 2) |
| y_acid | [0.4, 1.1] mol / mol | sec. 2 item 11 |
| Ea of the exits | [65, 85] kJ/mol, one value; the O2 exits' Ea drawn independently is a form | Deng 2021 |
| k_tr (trap) | [0.1, 10] /min | unsourced; numerical ends |
| t0 (fed-pot heat-up) | forms {0, 5 min} | unprinted |
| s_FOX | log-uniform within the measured DFA donor range (glucosone to MGO) | Hofmann Tables 3 - 4 (calibration rows, not scored) |
| 1-DG as a donor: s_DG1, phi_DG1, kappa_DG1, k0_DG1 | sec. 4.3 | sec. 4.3 |

### 4.2 The k_d band (sourced; label pH scale throughout)

All values are at 100 C, on the initial 25 C label pH. They are moved to label pH 7.0 by s_pH in [0.20, 0.54]:

| source | value at its own conditions | at label pH 7.0 (arithmetic) |
|---|---|---|
| Martins, fed DFG | 0.0368 /min at pH 6.8 | 0.0404 - 0.0472 |
| Deng, Met-ARP in 0.1 M phosphate | 0.021 /min at pH 7.5 | 0.0113 - 0.0167 |

- **From 0.1 M to 0.5 M phosphate.** The factor runs from 1 (saturated) to 2.8 (linear). The linear end uses Deng's
  water-to-0.1 M ratio under EX-058's non-saturating law. That law was measured on the amine + sugar condensation
  (Potman), so **for ARP decay it is a declared transfer to another step**, and it carries item 2's pH confound.
- **Band: [0.0113, 0.132] /min, log10 [-1.95, -0.88].** The two amino acids' spread stands in for the Phe transfer.
- **The floor form, 0.005 /min (log10 -2.30), is required agreement** (M7). The trapped pot's text end point implies
  ~0.006 - 0.008 at 100 C (sec. 5), below the band.
- **The band ends are transfer-class, not measured-class** (M7). A verdict whose profile depends on either end does not
  ship (sec. 11).
- **Starts for the mode searches.**
  - start 0: the centre, log10 -1.41;
  - start 1: Martins' and Deng's geometric mean, log10 -1.61;
  - start 2: the lower quarter, log10 -1.68.

### 4.3 The unmeasured donor's efficiency (1-DG; R2's main freedom)

- **s_DG1** is drawn log-uniform on [s_min, s_max]:
  - s_min is the lowest measured DFA share (glucosone, Table 4);
  - s_max is the highest (MGO, Table 3 plateau) x 10^0.82 (Granvogl; generous to R2, sec. 2 item 8).
- **Forms (required agreement for an R1-over-R2 verdict):**
  - s_max x 10^1.24 (CF-073's largest per-donor gap);
  - **Voigt's O2 factor**: s_DG1 x U[1.5, 1.6] in air pots (M2).
- **phi_DG1** is in [0.14, 0.30], the measured phi_D range.
- **kappa_DG1** is set as for the other unmeasured donors.
- **k0_DG1** is in [0.07, 1] /min at 100 C. The lower end transfers Davidek 2006a's 120 C bound; the upper end is
  numerical. k0_DG1 also sets 1-DG's trap efficiency (sec. 3.2).

### 4.4 The rivals

| id | what differs | new coordinates | band (source) | starts 0 / 1 / 2 |
|---|---|---|---|---|
| **R0-pinned** (the current model at shipped constants) | no ARP exits to donors, no route B, no acid. Route A on GO and MGO only, at B18's **shipped** k_go_ak / k_mgo_ak (log10 -6.54 / -7.53, read from `parameters_pyrazine.FROZEN_B18`), at the engine's shipped pH factor at label 7.0, with Phe charged as glycine | -- | -- | -- |
| **R1, the claim** | one k_d,7 and one share vector in every pot; g_B > 0; one psi | sigma_ox solved by 1-D root-finding so that the fed ARP's air excess total (PA + PAA) is matched; psi | sigma_ox in [0, 1); psi in [0, 1] | sec. 4.2; psi 0.5 / 0.35 / 0.65 |
| **R2, dicarbonyl only** (the case's R2, widened, M2) | g_B = 0. The fed ARP's air excess comes from eG + eFox and route A on their osones; its argon PA from 1-DG, 3-DG, fragments and released Phe | sigma_ox solved as for R1, with g_B = 0 | as R1; s_DG1 and its forms (sec. 4.3) | as R1 |
| **R3, pot-specific decay** | R1 with k_d,7,fed and k_d,7,glc drawn independently from the same band; shares shared. **R3 is identical to R1 on every fed-pot unit** | log10 k_d,7,glc | sec. 4.2 | crossed |
| **R4, the claim with a decay law** | R1 with order n (rate k_d (ARP / 100 mM)^(n - 1) ARP) and the in-situ pH term s_pH (sec. 3.4), **active in every pot** | n; s_pH; y_acid | n in [1, 2] (measured orders: Martins and Deng first, Wei 2022 second); s_pH and y_acid as sec. 4.1 | n 1 / 1.5 / 2 |

**Feasibility (MINOR, made non-circular).**
- For R1 and R2, sigma_ox is the root of "the predicted air excess total = the measured one", given every other
  coordinate of the draw. A draw with no root in [0, 1) has weight 0, and the share of such draws is reported.
- In R1, a drawn g_G + g_F that already over-produces the air excess (which would need g_B < 0) is infeasible, weight 0.
- c_Cu5 is inverted as in sec. 3.3, inside its bound.

**What the structure implies (M4).**
- **R4 contains R1.** At n = 1, R4 differs from R1 only through the pH drift, which is <= ~0.2 dex in the 0.1 M pots
  (sec. 3.4). **sigma-ARP cannot kill R4.** R4 loses only on Occam (its integrated n, s_pH and y_acid), or wins where
  order or drift matter: U3, then U-ARPglc.
- **In U3, order and pH drift are confounded** (M5). The 0.5 M pot raises k_d through n > 1 and lowers it through
  acidification.
- **R3 equals R1 on the fed-pot units.** So the route contrast R3 vs R2 is the R1 vs R2 result reused, and **outcome
  rows 2 and 3 reuse the R1 / R2 route verdict.** The power check has **three** pairs.
- **R0-pinned** represents only the glucose pot's PA. On every other unit it reads "cannot represent", a structural zero
  that is never summed.

**Nulls**, scored beside every rival:

| unit class | null |
|---|---|
| ratio units | constant ratio: the fed donors' mean under the same atmosphere |
| shape units | the typical shape: the fed MGO Table 3 shape |
| level units | typical yield: the fed pots' mean yield x the integrated forcing |
| **U3** | atmosphere-free: each species' yield equal in argon and air + Cu (the two arms' geometric mean, per time); Phe liberated linear in time through the arm mean |
| **U-ARPglc** | k_d,glc drawn from the band prior, with no transfer (R3's glucose marginal) |

## 5. Pre-screen arithmetic (disclosed; no band, start or tolerance comes from it)

1. **The fed ARP's early shape (D1).** The first-order k matching PA(10)/PA(30) = 0.381 is 0.014 /min; the one matching
   PA(500)/PA(120) = 1.61 is 0.0078. The early rise is faster than a single exponential. A product-of-two-releases R2
   predicts ~0.04, while free impurity Phe and a fast osone steady state allow ~0.1 - 0.4. Not dry-decidable.
2. **The glucose-pot bound.** k_d,glc >= 0.0126 - 0.0144 /min, weakened by sec. 6's noise forms.
3. **The trapped pot's Phe (dossier text, SEEN).** "~55 % at 320 min, 85 C" gives k_d(85) ~ 0.0025 and k_d(100)
   ~ 0.006 - 0.008 at Ea 65 - 85. **Acidification (sec. 3.4) can explain 0.07 - 0.7 dex of this**, and so can order:
   R4 / B1.
4. **R2's stoichiometry.** The air excess PA is 4.9 umol/mmol and the air + Cu excess 13.2. At DFA's measured donor
   efficiencies, R2 needs >= ~0.5 - 0.7 mol free osone per mol ARP in air and ~1.0 - 1.9 in air + Cu. Only an osone
   >= ~0.3 dex more efficient than any DFA fed donor rescues it. That is now FOX at most at MGO's efficiency, 1-DG up
   to x10^0.82, and the Voigt and CF-073 forms.
5. **O2 demand vs supply (M1; my arithmetic; a 1 mmol pot in a 10 mL closed vial).**

   | rival | O2 needed at 120 min, air | at air + Cu 5 mM | headspace air that supplies it |
   |---|---|---|---|
   | R1 (route B events only: air excess PA + PAA 7.6 umol; Cu 20.5) | 3.8 - 7.6 umol (+ any eG / eFox it draws) | 10 - 21 umol | **5 mL (43 umol) suffices** |
   | R2 (>= 0.5 - 0.7 mol osone per mol ARP; Cu 1.0 - 1.9) | **250 - 700 umol** | 500 - 1 900 umol | 29 - 82 mL in air; 58 - 222 mL in Cu |

   A closed vial holding 10 mL of solution rarely has 30 - 220 mL of headspace. **Under the 5 and 20 mL forms, R2's air
   excess is infeasible** at measured efficiencies. Under 100 mL (or pure O2) it is feasible. Hence these are
   required-agreement forms, never the primary.
6. **The acid load per pot**: sec. 3.4's table.
7. **The fed ARP's PA / PAA drifts 3.33 (10 min) -> 1.11 (500 min)** in Fig. 1 (arithmetic).
   - R1's route B has one psi, so this drift must come from a late, acid-rich route A: released Phe on ARP-made osones
     carrying the omega acid channel. Table 5 rules out PA -> PAA.
   - R2's drift follows its donor mix.
   - So U4's Fig. 1 ratios test R1's single-psi structure as much as R2. This is declared, and it shaped sec. 16.

## 6. Units: value basis, error terms and role (rules 1 - 3)

**Error model, common to every unit.**
- **Omega = 0.419 dex per predicted point, a declared transfer** (H0 EX-059, warm set; the larger of 0.419 / 0.390).
  - Its forms are 0.5x / 2x. The shared-reference covariance is omega^2 (I + 11^T); a ratio carries 2 omega^2. The
    within-pot rho is 0 (primary) or 0.5 (form).
  - **What omega does to power.** U4's largest expected difference (R2's glucosone fingerprint ~3.3 vs 1.8 measured,
    0.26 dex) is far below 2^0.5 omega = 0.59, so **U4 has essentially no power under omega**. It is kept as a claim
    unit and expected to tie. sigma-fp and every shape signature's screen tolerance (sec. 8) are vacuous; the screen's
    kills come only from the feasibility signatures, which carry no omega.
- **Detection floors (rule 1).**
  - PA / PAA: every decisive value is quantified, so there is no floor.
  - **U3: a keyed marker at or near zero is left-censored at the axis read resolution**, ~0.003 - 0.005 mmol/mmol ARP,
    the value set by the keying record (M10). The x2 form is required agreement. Deng's LOD is not used.
  - A floor never clips a prediction and is never set at a data point.
- **Between-laboratory tau (rule 2):** none in phase 1 (all within DFA). A cross-run level carries sigma_run =
  0.075 dex (the 7 - 18 % run scatter of CF-074 / 075).
- **Vessel per time point (M6).** One vessel per sampled time cannot be excluded for any Hofmann series. **Form "sigma_run
  between points"**: an independent 0.075 dex on every point of U1e, U1l, U5a, U5b, U-ARPglc and U3 (in addition to
  each run's shared offset). Required agreement.

| unit | data | value basis | measurement noise | other terms | floor | role |
|---|---|---|---|---|---|---|
| **U1e: fed early shape** (D1) | Fig. 1 circles: log10 PA(t)/PA(120), t = 10, 30, 60 | within-run ratio | read error; rounding +-0.05 umol/mmol (sd 0.029); s_rep 0.025 | omega (shared reference); the between-points form | none | claim, route |
| **U1l: fed late shape** | the same, t = 240, 300, 500 | as U1e | as U1e | as U1e | none | guard |
| **U4: fingerprints** (D4) | log10 PA/PAA: Table 2 fed ARP (argon, air, air + Cu); Fig. 1 circles 10 - 500 min | within-vial ratio | print +-0.05 per analyte; s_rep; Fig. 1 read + rounding | 2 omega^2 | none | claim, route (no power under omega: expected tie) |
| **U3: trapped osones and Phe** (D3; yields blind, Phe end point and ordering seen) | Figs. 2 - 3 keyed: 1-DG, 3-DG, glucosone and Phe, 20 - 320 min, argon and air + Cu 25 mM, 85 C | **absolute yield per mol ARP charged** | keyed read error; calibration coordinates N(0, 0.10) (3-DG, glucosone), N(0, 0.20) (1-DG; form 0.30), N(0, 0.05) (Phe) | omega; **sigma_run, one offset per arm** (each arm is one vessel run, its own run vs Table 2's calibration vials); the arp0 band; the between-points form | censored at the axis read resolution | **claim, route (primary); claim, decay law** (M5, after B1's pH form) |
| U3a, U3b (conditional, sec. 12.1) | U3a: the argon arm's trapped 1-DG yield against the argon fed PA (Table 2, 100 C); U3b: glucosone trapped per Phe liberated within Fig. 3 | U3a cross-pot ratio (sigma_run); U3b within-run ratio | as U3 | 2 omega^2 | as U3 | added before the freeze **only if** U3's synthetic power fails; both free of the Cu law |
| **U5a: glucose pot, hofmann2000 Fig. 2** | PA and PAA, 10 - 500 min, as log levels | absolute within lab, forced by another run | read errors; s_rep | sigma_run; omega; the between-points form | none | claim, transfer and decay law |
| U5a-PA-shape | log10 PA(t)/PA(120) from the same points | within-run shape | as U5a | omega | none | **R0-pinned's shape-only comparison** (a favoured rival must beat R0-pinned here too) |
| **U5b: Fig. 1 squares** | as U5a (the 10-min squares kept with their 0.126 / 0.127 dex error) | as U5a | as U5a | as U5a | none | claim, transfer and decay law; its PA shape also as an R0-pinned shape unit |
| **U-ARPglc** | Table 3 ARP-Phe: each interval's ln(ARP(t1)/ARP(t2)) <= k_d,glc,eff (t2 - t1) | censored, scale-free | integer rounding (uniform +-0.5 umol) **plus the analyser CV integrated over [2, 10] %**; form: sigma_run between points | **form "isomers"**: the measured ARP = ARP + I, with I inert, formed from ARP and reaching f_iso in [0, 0.3] of the measured value at 300 min | none | claim, transfer and decay law |
| G-T3 | hofmann2000 Table 3 shape (T1, T1b) | within-run | as EX-049 | -- | none | precondition and Holm guard |
| G-glcO2 | Table 2 glucose pot PAA, argon and Cu | -- | -- | -- | -- | reported (not forceable; the PA rows were B29 fit rows) |
| R-D2X, R-Deng | Martins fed vs in-pot DFG ratio; Deng +-Met | ratio | HPD; read error | transfer | ARP LOD 0.18 mM (Deng) | reported only (band sources) |

**Conditioning (prequential; both blocks ordered; declared now).**
- **The fed block** is scored in the order U4 -> U1e -> U1l -> U3. Each unit is conditioned on the calibrations and the
  units before it.
- **The glucose block** (U5a, U5b, U-ARPglc, in that order) is conditioned on the calibrations plus the fed block, in
  two forms:
  - **(F) full:** all four fed units;
  - **(T) "100 C fed only":** U4, U1e, U1l and Table 2, with **U3 excluded**.
- The transfer and decay-law contrasts require agreement between (F) and (T). If they disagree, the reading is
  "confounded by acidification" (B1).
- **The total over a block does not depend on the order** (the chain rule). The per-unit attribution does, and the
  order above fixes it.

## 7. Keying hofmann2000b Figs. 2 - 3, blind

- **Done.** Two independent agents digitised every plotted series in Fig. 2 (argon) and Fig. 3 (air + CuSO4 25 mM),
  as briefed by the coordinator. The series are the 1-DG, 3-DG and glucosone quinoxalines, Phe liberated, and ARP-Phe if
  plotted, at each plotted time, with the y-axis quantity, unit, caption and Methods verbatim. Markers near zero are
  flagged, never omitted. Their files sit in the coordinator's scratchpad under `hofmann2000b_keying/`. **The author has
  not opened them.**
- **Adjudication: a third agent with no design access** (M10) applies the established convention: mean(A, B); read
  error max(err_A, err_B, |A - B| / 2); exposure recorded in the file header. That agent also records the axis read
  resolution, which sets U3's censoring level.
- **Commit order.** The adjudicated `data/lit/fit_rows/strecker_hofmann2000b_fig2_fig3.yml` (with the reader-B copy and
  `source_reading_reconciliation_strecker_hofmann2000b_fig2_fig3.md`) is committed **only after the v2 freeze commit,
  on top of it**. The ledger update is the register session's (one writer).
- **V-blind (enforced in the generator)**, before any real-data trapped-pot run:
  - (i) the frozen v2 commit (sec. 20) is an ancestor of the keyed-file commit (`git merge-base --is-ancestor`);
  - (ii) the keyed-file commit is an ancestor of HEAD;
  - (iii) the keyed file's hash at scoring equals its hash at its commit;
  - (iv) no trapped-pot prediction on real data is printed or stored before (i) - (iii) hold.
- **Before then,** the trapped-pot sub-model runs only on synthetic inputs (V3, sec. 12.1).

## 8. Mechanism screen (rule 9; before the freeze; >= 600 Sobol draws per rival)

- Every rival, every integrated coordinate drawn, the calibrations inverted per draw. Coverage is reported when fewer
  than 10 draws pass.
- **The screen runs on every signature except sigma-trap, before the keyed file is committed.** sigma-trap runs after
  the keying commit, as part of scoring, and can only demote.
- **Tolerances.**
  - Feasibility signatures are hard constraints.
  - Shape signatures use z_0.975 sqrt(measurement variance + omega^2 c). **Under omega these are near-vacuous**
    (sec. 6), so a shape signature cannot kill in practice. That is declared.

| sig | data | measured | bears on |
|---|---|---|---|
| sigma-stoich | Table 2 fed ARP, air and air + Cu: sigma_ox has a root in [0, 1) and c_Cu5 lies in [1, 25] (sec. 4.4) | -- | R2 above all |
| sigma-O2 | the same under each finite-budget form (sec. 3.3) | -- | R2 (reported per form; a kill needs every form, sec. 8 kill rule) |
| sigma-ARP | Table 3 intervals vs k_d,glc,eff under the noise forms of sec. 6 | 0.0126 - 0.0144 at 30 -> 60 min | R1 (R3 meets it by construction; **R4 cannot be killed by it**, since it contains R1) |
| sigma-lag | Fig. 1: log10 PA(10)/PA(30) | -0.42 | R2 (vacuous under omega) |
| sigma-fp | Table 2 fed ARP PA / PAA in air + Cu | 1.8 | R2 (vacuous under omega) |
| sigma-cross | glucose pot: PAA overtakes PA between 60 and 120 min | yes | all |
| sigma-early | hofmann2000 Fig. 2 at 30 min | -- | all |
| sigma-trap | Fig. 3 (keyed): trapped glucosone within tolerance at the sigma_ox the calibration needs | after keying | R2 |
| sigma-R0 | R0-pinned's glucose-pot PA | printed | -- |

**Kill rule.**
- A rival is killed only if no draw reproduces every signature that bears on it, confirmed by EX-049's bounded violation
  search, under **every** form that bears on it. For R2 those are the CF-073 and Voigt forms, o_Ar = 0.1, the 100 mL
  headspace and the unbuffered trapped pot.
- **A screen kill of b credits a only if a passes the adequacy gate (sec. 9.4) on that contrast's claim units** (M3).
  Otherwise the outcome is row 5.

## 9. Densities, tie band, required agreement, adequacy, multiplicity (rules 4 - 8, 11)

### 9.1 Predictive density per rival, then differenced

- **The density.** lppd_uk is the log of the integral of the unit's likelihood over the band prior, conditioned per
  sec. 6, computed separately for each rival and then differenced.
- **The estimator.** PSIS with a unit-tuned Student-t proposal (5 df) at the mode, found by <= 300 evaluations from the
  three starts. Accepted iff k-hat < 0.7, ESS >= 100 and SE_MC <= 0.15 sd_rep / sqrt(2). S runs 256 -> 4 096, then
  Gauss-Hermite on the coordinates that carry > 90 % of the linearised variance (budgeted, sec. 15). Failing that, the
  unit reads "density not resolved" (a tie).
- **Every form gets an exact density.**
  - **Likelihood-only forms reweight the same draws**, with k-hat and ESS re-checked: omega 0.5x / 2x, rho, the
    censoring x2, the analyser-CV band ends, the 1-DG calibration sd, sigma_run between points.
  - **Layer-changing forms get new draws:** o_Ar; CF-070; early-step forcing; t0; the Cu law; independent O2-exit Ea;
    CF-073; Voigt; the headspace x nu_O2 forms; pH drift on (trapped pot, non-R4 rivals); the unbuffered trapped pot;
    the k_d floor form; the isomers form; per-pot O2_reflux (R3); the +0.5 kappa ceiling.
- **Diagnostics.** Plug-in and linearised densities are printed, and any gap larger than the tie band is flagged.
- **V-density** on a two-coordinate R1 / R3 toy, before any score.

### 9.2 Tie band

- sd_rep(d_u) comes from 4 000 redraws of y_u from measurement noise only: read error, rounding, s_rep, sigma_run,
  calibration coordinates, and the analyser CV at its band midpoint. Omega is never part of it.
- The band is b_u = z_0.95 sd_rep.

### 9.3 Required agreement

A unit favours a iff d_u > b_u under the primary density and keeps its sign under every form of secs. 3 - 6 that applies
to it:
- **unsourced choices:** o_Ar, t0, the Cu law, headspace x nu_O2, the 1-DG calibration sd, the k_d floor, the isomers;
- **declared transfers:** omega, rho, censoring, CF-073, Voigt, independent Ea, CF-070, early-step forcing, the kappa
  ceiling, pH drift on, the unbuffered trapped pot;
- **for the transfer and decay-law contrasts:** both conditioning forms (F) and (T);
- **for R3 verdicts:** a per-pot O2.

**The same forms apply to "adequate"** (MINOR, as M1b): a rival counts as adequate on a unit only if it is adequate
under every form.

**A favoured rival must also beat R0-pinned** on U5a-PA-shape and U5b's PA shape, and **beat all three nulls** on each
claim unit.

### 9.4 Adequacy gate

- **The test.** A unit may favour a only if a is adequate there: posterior-predictive p > 0.05 for the Mahalanobis
  discrepancy under the full Sigma_u, or the censored PPC.
- **The alternative.** Failing that, the predictive means must differ by more than omega sqrt(c).
- **Otherwise "both inadequate".**
- **The record states, per favouring unit, which of the two conditions it passed.** It also states whether the verdict
  would change if the |Delta mu| alternative were removed (MINOR).

### 9.5 Preconditions and Holm

- **Preconditions.** Each surviving rival has a feasible calibration in >= 10 % of its draws. G-T3 passes in-sample. No
  claim unit is left "not resolved", or that unit is a tie.
- **Guard groups.** U1l and G-T3 for every pair; U1e and U4 also for transfer and decay law.
- **Holm's step-down at family-wise 0.05, one-sided**, across a pair's guard groups.
- **A lost guard** must name a candidate compensating error: the case's sec. 1 lists in phase 1, sec. 14.4 in phase 2.
  It must not outrank the best claim unit won; otherwise "credited, held".

## 10. Contrasts, decision rule, outcomes, stop rule, hour cap

### 10.1 Contrasts (three pairs)

| contrast | pair | claim units | guards |
|---|---|---|---|
| **route** | R1 vs R2 (**reused for R3 vs R2**: R3 = R1 on every fed unit) | U3 (primary), U1e, U4 (and U3a / U3b if added) | U1l, G-T3 |
| **transfer** | R1 vs R3 | U5a, U5b, U-ARPglc, under (F) and (T) | U1e, U1l, U4, G-T3 |
| **decay law** | R4 vs R1 | U3 (pH-drift form active in both), U5a, U5b, U-ARPglc | U1e, U1l, U4, G-T3 |

**Rule.** A contrast credits a over b iff all of the following hold:
- the primary claim unit (for transfer and decay law: the summed block with its own tie band) favours a;
- no claim unit favours b;
- Holm's guard step holds;
- adequacy and required agreement hold on each favouring unit.

A screen kill credits only under sec. 8's M3 condition.

### 10.2 Bound-active coordinates

Every verdict also passes sec. 11. A verdict conditional on a numerical or transfer-class bound does not ship.

### 10.3 Outcomes: exclusive, decided in this order

| order | condition | verdict | ships |
|---|---|---|---|
| 1 | route credits R1 over R2; transfer credits R1 over R3; decay law does not credit R4 over R1; R1 beats R0-pinned and the nulls | **"claim supported (single laboratory, PROVISIONAL)"** | the phase-1 layer is validated for DFA's design class, **relative within-lab only**. No engine change. Phase 2 unlocked |
| 1b | the route condition holds and decay law credits R4 over R1 (R4 over R2 follows from R1's route result where U3 does not reverse it under R4's pH law; checked) | "two routes; the decay moves with order / pH" | as 1, with R4's law |
| 2 | **the route result of row 1** and transfer credits R3 over R1 | "two routes; the fed decay does not transfer" | nothing for a new pot; phase 2 unlocked with the network's own in-pot ARP decay |
| 3 | **the route result of row 1** and transfer reads tie, both inadequate, "confounded by acidification", or demoted | "two routes; transfer undecided" | as 2 |
| 4 | route credits R2 over R1 (and hence over R3) | "dicarbonyl only" | nothing; the W6 candidate (sec. 19) |
| 5 | anything else | **the stop rule: "not separated on disk"** | nothing |

Rows 1 - 4 exclude one another through their route and transfer conditions, and row 5 is the complement.

### 10.4 Stop rule and caps

- **Stop rule.** After the screen, the power check and one scoring pass: if no row 1 - 4 condition holds, row 5 applies
  and phase 1 ends.
- **Author-work cap: ~22 h** (build 9 - 11 h, V-checks and the PSIS pilot 3 h, screen and power checks 2 - 3 h,
  scoring and record 3 h, review answers 3 - 4 h).
- **Compute cap: 24 h wall** (sec. 15).
- At either cap, row 5 applies unless a row 1 - 4 condition already holds.

## 11. Compensating structure: profiles (rule 12)

**What is profiled.** Every bound-active coordinate, at 0.25 / 0.5 / 1.0 dex (logs), 0.05 / 0.1 / 0.2 (fractions) and
1 / 3 / 10 kJ/mol inside the bound, with each decisive d_u recomputed. These include:
- O2_reflux at 1;
- k_d,7 at either end;
- cal_ARP;
- s_DG1 at s_max;
- k0_DG1 at 1;
- k_tr at either end;
- c_Cu5 at 25;
- n;
- s_pH and y_acid at their ends;
- sigma_ox near 1.

**Switched-off coordinates** get one extra point: g_B = 0 in R2 at R1's median; R3's k_d,glc at R1's k_d.

**What ships.** A verdict that depends on a **numerical** bound (k0_DG1 = 1, k_tr, c_Cu5 = 25, the kappa ceiling) or a
**transfer-class** bound (k_d,7's ends, s_DG1's upper end, s_pH's and y_acid's ends) does not ship (M7). Definitional
bounds (fractions, orders 1 / 2) do.

## 12. Power checks, then freeze (rule 10)

### 12.1 U3 first, on a synthetic design (before the keyed file is committed)

- **The design.** Both arms at 20 / 40 / 80 / 160 / 320 min; 1-DG, 3-DG, glucosone and Phe; read errors at the Fig. 1
  keying's level (<= 0.011 umol/mmol, plus 0.03 dex) and censoring at 0.004 mmol/mmol.
- **Truths.** R1 and R2 draws taken from the band prior **reweighted by the non-U3 fed units only** (U4, U1e, U1l,
  Table 2). Never the screen's best draw on U3 (M8).
- **The check.** 200 synthetic data sets per truth, scored by the route rule on U3 alone, under the primary and the Cu-law
  forms.
- **If U3 cannot separate,** U3a and U3b (sec. 6) are added before the freeze, with their own synthetic check. **Not
  run for v2**, because the layer is not built (sec. 21).

### 12.2 All three pairs

- **Truths.** Prior draws, or draws reweighted by units outside the pair's claim set (never a pair's own claim units).
- **The check.** R = 200 synthetic sets from each unit's full Sigma_u.
  - **U-ARPglc's generative model (M8).** For the truth's k_d,glc, take the formation F(t) >= 0 of minimum norm that
    reproduces the measured interpolant under that k_d, clipped at 0. Then ARP(t) solves dARP/dt = F - k ARP, and is
    rounded to integers with the analyser CV drawn from its band.
  - U3 is generated at the keyed times only after the keying commit; before it, at the synthetic design.
- **Scoring the synthetic sets.** Both rivals, by reweighting. **Precomputed draw sets** are kept for the primary and the
  three forms most likely to flip a sign (o_Ar, the Cu law, headspace 20 mL). Power on the other forms is labelled
  "primary-only, optimistic".
- **The binomial demotion rule (M8).** A pair is demoted unless power exceeds the false-credit rate by a one-sided exact
  two-proportion test at 0.05 (R = 200 each). A power below 0.5 is labelled "low power".
- **Then phase 1 freezes**, with the commit hash entered in sec. 20.

## 13. Order of operations

1. Re-review of v2 if the coordinator asks; then answers.
2. Build `m3_strecker_layer.py` and its tests:
   - **V1:** EX-049 bit-exact at its structure limit.
   - **V2:** mass and O2 balance to 1e-9.
   - **V3:** trapped-pot totals on synthetic inputs.
   - **V4:** R0-pinned from the code.
   - **V5:** the pH scale printed on every transfer.
   - **V-acid:** sec. 3.4's table reproduced by the layer's pH routine within 0.02 pH.
   - **V-time; a one-unit PSIS pilot at full dimension** (U3 under R2, the largest coordinate set): k-hat, ESS, wall
     time. The compute plan is re-estimated if it exceeds the cap.
   - **V-blind; V-density.**
3. The U3 synthetic power check (12.1), then the screen without sigma-trap (8), then the power check (12.2), then
   demotions only. **Freeze v2 (commit).**
4. The coordinator commits the keyed Figs. 2 - 3 on top. V-blind (i) - (iii) pass. sigma-trap runs (demotion only).
5. Scoring: every rival, every form, then profiles. No decisive number is printed before every density exists.
6. Results review by an independent agent; the outcome row; the EX entry; the phase-2 amendment if unlocked.
7. After any screen, only demotions are allowed.

## 14. Phase 2 (pre-declared conditional): the layer fed from the network

### 14.1 When and where

- Phase 2 runs after outcome row 1, 1b, 2 or 3.
- It edits `network.py` and the trunk lane, so under the one-writer rule it runs before A2 + H3 ports the trunk or on
  the ported module (the coordinator's scheduling).
- Amendment A-P2 freezes its bands and units before any run, with its own screen and power check.

### 14.2 Structure (declared now)

- **New species:** PHE, PHEARP, PHEAL, PHEAC (provisional keys); DG1 and FOX if the trunk lacks them. Leu and Met are
  identity-scaled and reported only.
- **New steps:**
  - Glc + Phe -> PHEARP: the trunk's glycine Schiff / Amadori constants with an identity term (A-P2);
  - the PHEARP exits of sec. 3.2: phase 1's posterior (rows 1 / 1b), or the network's own in-pot decay (rows 2 / 3),
    with R4's law if row 1b;
  - per-amino-acid Strecker steps at phase 1's constants.
- **B18's `r_go_ak` / `r_mgo_ak` keep their constants.** The aminoketone co-product is booked to AKG / AKM (primary) or
  to FRAG_N (form).
- **The non-Strecker amine sink (CF-095) is M6's** and is not added.
- **pH** is the network's balance or the MTO, with the scale stated.

### 14.3 Phase-2 units (all scale-free)

| unit | data | basis | tau | role |
|---|---|---|---|---|
| U2-Hof | the glucose-pot shapes and Table 3 ARP / Phe, from the network's own courses | within-run | none | claim |
| U2-Weenen | Weenen 2001 Table II pH 7.0, x4.2 | cross-lab within-study ratio | none | claim |
| U2-Deng | Deng 2022 methional ARP / sugar ratio, 100 and 120 C | ratio | none | reported (engine-seen) |
| -- | any absolute cross-lab level | absolute | tau_ABS 0.964 + figure term | reported only |

### 14.4 Side effects (rule 15)

**The owner signs off at ship time.** Before and after, every row is re-scored: the full held-out contrast benchmark
(coverage first); B18's 10 fit rows; the trunk's dicarbonyl rows; the W1.1 methionine and Strecker items. Cremer 2000,
Balagiannis 2009, Huang 2017 (TUM Becker), Meynier 1995 Table 4 and the case's other TEST cells are read as item 7 only.

| expected movement | why | candidate compensating error |
|---|---|---|
| **every pyrazine row: up** | the new Strecker steps make AKG / AKM (primary form); under FRAG_N, only donor depletion acts | B18's pyrazine-yield constant (EX-021; M6's 2.9-dex deficit) |
| trunk MGO / GO / 3-DG / glucosone in Phe / Leu / Met pots: down | consumption by the Strecker steps | the trunk's MGO level (B21); the missing 3-DG / glucosone donors (F16-b) |
| methionine items and B22's rows: from refused to predicted | methional gets a route | B22's refused constants; glycine as the only amine |
| Leucine requests: answered | coverage | -- |
| pH in Phe pots: slightly more acid | PAA charged as a carboxyl | the B51 acid yield (bound-active at B60) |
| sulfur and acrylamide lanes | only if they read the trunk's dicarbonyls: **UNVERIFIED, grepped at the A-P2 build** | -- |
| thiol rows | none expected | a loss is a build error until shown otherwise |

## 15. Compute (wall-clock hours, <= 8 cores, `nice`d, memory checked at launch, never paused)

**Phase 1 runs on the standalone layer.** Assumed **50 - 150 ms per draw** serial (M9), covering ~6 ODE pots plus O2 and
pH states; UNVERIFIED until V-time. The machine delivers ~3.3 serial-equivalents at 8 workers (profile 2026-11-05).

| phase 1 step | count | serial | wall |
|---|---|---|---|
| U3 synthetic power (12.1) | 2 truths x 200 sets x reweighting of ~4 096 draws x 2 forms (draws ~16 k) | 0.2 - 0.7 h | 0.1 - 0.2 h |
| screen | 5 rivals x 600 - 1 200 x 0.05 - 0.15 s | 0.04 - 0.25 h | ~0.1 h |
| mode searches | ~7 units x 5 rivals x ~11 forms x <= 300 | 1.6 - 4.8 h | 0.5 - 1.5 h |
| **PSIS** | **~7 units x ~11 forms x 5 rivals (385 unit-form-rival cells) x S 1 024 - 4 096 x 0.05 - 0.15 s** | **5.5 - 66 h** | **1.7 - 20 h** |
| Gauss-Hermite fallback | assume 20 % of cells, 3 coordinates x 10 nodes (1 000 evaluations) | 1.1 - 3.2 h | 0.3 - 1 h |
| power check (12.2) | 3 pairs x 2 truths x 200; precomputed draw sets for 4 forms | ~1 - 2 h | 0.3 - 0.6 h |
| profiles | ~14 coordinates x 3 points x <= 1 200 | 0.7 - 2.1 h | 0.2 - 0.6 h |
| **phase 1 compute** | | **~10 - 80 h** | **~3 - 20 h wall (credible range)** |

- **Compute cap: 24 h wall.**
- **Declared fallback** if the PSIS pilot projects more than the cap:
  - (i) S is capped at 1 024;
  - (ii) forms whose reweighting fails k-hat are run at S = 1 024 only;
  - (iii) any cell still unresolved reads "density not resolved" (a tie on that unit);
  - (iv) the record says "compute-capped".
- **No form is dropped to fit the cap.**

**Phase 2 runs on the current engine** (profile: median pot 0.068 s, slowest 0.46 s; ~16 pots per draw at ~0.15 s, so
~2.4 s per draw; UNVERIFIED until A-P2). Screen ~0.4 h wall; PSIS 3 rivals x S 1 024 - 4 096 x 3 forms ~1.9 - 7.5 h;
power check ~0.3 h; side effects 1 - 4 h (UNVERIFIED scorer time). **Total ~3.6 - 12 h wall**, on top of 24 - 30 h of
author work.

## 16. Expectations (written now; not blind; coherent across rows)

| event | P |
|---|---|
| sigma-stoich or sigma-O2 kills R2 under the primary bands / under every form (incl. 100 mL, CF-073, Voigt) | 0.35 / 0.08 |
| the screen kills R4 | 0 (it contains R1) |
| R1 survives the screen / R3 survives | 0.75 / 0.9 |
| U3 synthetic power >= 0.5 under the primary / under the Cu-law forms too | 0.6 / 0.4 |
| route: R1 over R2 / R2 over R1 / not separated | 0.40 / 0.05 / 0.55 |
| transfer (given the route result): R1 over R3 / R3 over R1 / tie, inadequate, confounded | 0.2 / 0.2 / 0.6 |
| decay law: R4 over R1 | 0.1 |
| outcome rows 1 / 1b / 2 / 3 / 4 / 5 | 0.07 / 0.03 / 0.08 / 0.22 / 0.05 / 0.55 |
| a pair demoted by the power check | 0.45 |
| the keyed Fig. 3 glucosone yield < 0.5 mol/mol at 320 min | 0.65 |
| phase 1 verdict by 12-05 / phase 2 shipped by 12-26 | 0.4 / 0.12 |

**The rows are coherent.** Rows 1 + 1b + 2 + 3 = 0.40 = P(route credits R1). Row 4 = P(R2 over R1). Row 5 = P(route not
separated) = 0.55. Row 1 = 0.40 x 0.2 x 0.9; row 2 = 0.40 x 0.2.

## 17. Owner-level decisions

1. **Shipping phase 2 moves shipped pyrazines and the trunk's dicarbonyl pools**: owner sign-off at ship time. Recorded;
   nothing to decide now.

**Closed by the coordinator (2026-11-05):**
- phase 1 ships relative within-lab only (one laboratory);
- no data request now: the fed 1-DG + Phe pot and the fed ARP-Phe pot with the ARP measured per atmosphere are W6
  candidates (sec. 19).

**Decided here, under the rules:**
- O2_reflux shared, with a per-pot form for R3;
- R3's two decays on one band;
- psi integrated;
- Martins and Deng as band sources, so their ratios are reported only;
- both blocks prequential, plus the (T) form;
- the Table 2 glucose PA rows diagnostic;
- the 10-min squares kept with their error;
- phase 2's bookkeeping forms;
- B18's glycine steps unchanged;
- every unsourced choice is a required-agreement form.

## 18. Adversarial review, and open questions

**Review.**
- (a) v1's design review is done (MORE FIXES), and v2 answers it (sec. 21). A re-review of v2 runs if the coordinator
  asks. It should check B1's buffer arithmetic, the finite-O2 forms, R2's widening, and whether omega makes the route
  contrast powerless; the power check should then demote it.
- (b) A results review follows scoring, before the outcome row is written.

**Open questions.**
1. What do Figs. 2 - 3 plot (units, ARP-Phe remaining, the 1-DG standard)? The adjudicated file answers this; it is
   read only after the freeze.
2. Does 0.6 M 1,2-diaminobenzene react with ARP-Phe or Phe, and does trapping change the acid load? It is modelled
   (sec. 3.4) but unmeasured.
3. Were the Fig. 1 circles and Table 2's fed rows made from one ARP-Phe batch?
4. Does Deng's Ea hold for Phe-ARP and for the O2 exits?
5. How does Cu act from 5 to 25 mM?
6. Does Table 3's "ARP-Phe" include isomers? This is now a form (sec. 6).
7. Phase 2's timing against A2 + H3's one-writer rule.
8. The charge and stripping in Weenen's pH 7 pots.
9. Do the hofmann2000 and 2000b reflux pots share one O2 regime?
10. Does H2's bake-off run on this layer after the freeze? Recommended: yes.
11. **What is Table 2's vial headspace?** It is unprinted, and it decides sigma-O2's forms.
12. **What was the trapped pot's phosphate: 0.5 mol/L or the Methods' "0.5 mmol/L"?** The unbuffered form covers the
    latter.

## 19. For the register session (nothing here is edited by this design)

**EX entry (id to be assigned; EX-061 is the highest on origin/main at cc5cbb55):** "M3 phase-1 pre-registration
(design, v2 draft)".

**CF candidates:**
- (a) **(downgraded)** At 0.5 M phosphate under a linear law (a transfer from EX-058's step), EX-049's in-pot 0.058 -
  0.069 /min is *consistent with* the measured amino-acid ARP decays: Deng ~0.057, Martins ~0.10. Deng's 1.8x phosphate
  factor carries a pH confound.
- (b) The glucose pot's ARP-Phe fall (30 -> 60 min) needs k_d >= 0.0126 - 0.0144 /min at face value. It weakens under
  analyser CV and vessel-per-time noise.
- (c) Fed decays look slow in DFA's fed designs: Fig. 1's late shape ~0.008; the trapped pot's text ~0.006 - 0.008. The
  trapped pot's acidification (0.07 - 0.7 dex) and the order are confounders.
- (d) EX-049's anoxic-route bookkeeping (note on EX-049).
- (e) R2's stoichiometric need (air ~0.5 - 0.7, Cu ~1.0 - 1.9 mol osone / mol ARP).
- (f) **NEW: R2's O2 demand** of 250 - 700 umol in air for a 1 mmol pot (29 - 82 mL of headspace air), against R1's
  4 - 8 umol. Table 2's headspace is unprinted.
- (g) **NEW: the acid load of the trapped pot** (0.5 M ARP in 0.5 M phosphate): a pH drop of 0.35 to >= 1.1, or an
  unbuffered pot under the Methods' "0.5 mmol/L". Its decay is not the 100 C pots' decay.
- (h) The omega transfer from thiols to Strecker species.
- (i) The Weenen x4.2 / Hofmann x3.9 agreement (a check, not a conflict).
- (j) **NEW: Martins 2003's "pH 5.5" arm drifted 0.8 - 0.9 units**, so its initial-label slope (0.69) overstates the
  drift-corrected one (0.54 - 0.56). This affects every use of Martins' pH slopes (M1b's bands too).

**Wording errors:**
1. Plan row M3's D2 sentence (as (a)).
2. Case sec. 2.4's "0.0062 floor" (as (b)).
3. EX-049's Table 2 air ratio was in-sample.
4. Case D6 is phase 2 only.
5. v1's sec. 2 item 3 ("D2 cannot reject R3" stated as a bound); v2 downgrades it.

**W6 candidates (coordinator: not requested now):**
- a fed 1-DG + Phe pot;
- a fed ARP-Phe pot with the ARP measured under argon and under air, with the headspace stated.

**Keying:** the adjudication (a third agent) and the commit of Figs. 2 - 3 after the v2 freeze; a
FIT_HOLDOUT_DECLARATION amendment for M3's rows.

**Tags:**
- Martins 2003 / 2003b and Deng 2021: the k_d and s_pH bands and y_acid;
- Davidek 2006a / 2005: y_acid and k0_DG1;
- Granvogl 2006: s_rep and s_DG1's upper end;
- Voigt 2010: the s_DG1 form.

**Results README:** a `DESCRIPTIONS` entry for this file in `scripts/generators/build_results_readme.py`, then the
README regenerated, before merge (`tests/unit/test_results_readme.py` fails until then).

## 20. Results

*(the freeze commit is entered here first; results appended after the run; nothing above this line is edited after
the first number exists)*

## 21. Revision log

- v1 DRAFT, 2026-11-05 (a66fe933): first version.
- v2 DRAFT, 2026-11-05, after the independent design review (MORE FIXES, no redesign) and the coordinator's decisions,
  before any build:

| finding | change |
|---|---|
| B1 acidification | per-pot acid-load table (my arithmetic, y_acid 0.4 - 1.1 from Martins 2003b and Davidek 2006a; sec. 3.4); R4's pH term active in phase 1; "pH drift on" a required form for the trapped pot; the (T) conditioning form (100 C fed units only, no U3) with required agreement against (F), else "confounded by acidification" (secs. 3.4, 6, 9.3, 10.3) |
| M1 O2 budget | finite-headspace forms {5, 20, 100 mL} x nu_O2 {0.5, 1}; the demand-vs-supply table per rival (sec. 5 item 5); sigma-O2 (secs. 3.3, 8) |
| M2 R2 widened | the oxidative fragment exit eFox (efficiency within the measured donor range), Voigt's x1.5 - 1.6 form on s_DG1, species-specific trap efficiency from competition with k0_DG1 (secs. 3.2, 4.3) |
| M3 | a screen kill credits only if the survivor is adequate on the contrast's claim units, else row 5 (secs. 8, 10.1) |
| M4 | R4 contains R1 (sigma-ARP cannot kill it; P(kill) 0); R3 = R1 on the fed units, so three power pairs; rows 2 - 3 reuse the route result (secs. 4.4, 10, 12, 16) |
| M5 | U3 a claim unit of the decay-law contrast with the pH form; order and drift confounded there (secs. 4.4, 6, 10.1) |
| M6 | U-ARPglc: analyser CV integrated over [2, 10] %; sigma_run between points (also on U1, U5, U3); the isomers form (sec. 6) |
| M7 | k_d band re-derived with Deng's 100 C slope (0.19 - 0.20) and Martins' drift-corrected slope (0.54 - 0.56): [0.0113, 0.132]; floor form 0.005; ends transfer-class; CF (a) downgraded; the 1.8x pH confound; EX-058's transfer declared (secs. 2, 4.2, 11, 19) |
| M8 | truths never from U3's best draw; U-ARPglc's generative model; precomputed draw sets for four forms, the rest labelled optimistic; a binomial demotion test; the U3 synthetic power check before the keyed commit, with U3a / U3b (free of the Cu law) as the pre-declared fallback (secs. 6, 12) |
| M9 | compute per unit x form x rival at 50 - 150 ms; Gauss-Hermite budgeted; a 24 h cap with a declared fallback; V-time plus a one-unit PSIS pilot before the freeze; credible phase-1 wall 3 - 20 h (secs. 13, 15) |
| M10 | the U3-Phe end point and the ordering marked seen; V-blind (i) - (iv) with the freeze commit an ancestor of the keyed commit and the hash unchanged; a third-agent adjudicator; censoring at the axis read resolution (secs. 0, 6, 7) |
| MINORs | both blocks prequential with the order stated; omega as a declared transfer, with U4, sigma-fp and the shape screen powerless under it; required agreement on "adequate"; the |Delta mu| route reported; coherent probabilities; row 1b in sec. 1.2; non-circular feasibility, wB < 0 infeasible, c_Cu bounded; U3 sigma_run; nulls for U3 and U-ARPglc; Martins' k1 / k2 / k3 labelled; Granvogl generous to R2; the R0-pinned shape-only units; Fig. 1's PA / PAA drift declared (sec. 5 item 7) |
| owner items | (2) closed: relative within-lab only; (3) closed: W6 candidates, no request now |

**Not applied, with reasons:**
- **The U3 synthetic power check was not run.** It needs the layer, which does not exist; this file is design only, as
  briefed. It runs at build, before the freeze and before the keyed commit (sec. 12.1).
- **The acid-load numbers were not computed by the layer.** They come from a standalone Henderson-Hasselbalch script,
  not the layer's own pH routine; V-acid checks the two against each other at build.
