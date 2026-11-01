# E6 + E11 pre-registration: MFT's two intact-C5 routes, refitted with per-route formation barriers

*Written 2026-11-01 (night), before any refit, member or stage-0 score existed. PROPOSAL for the owner's morning
decision. Nothing ships, nothing in `src/` is edited, nothing is committed. Origin:
`docs/reference/pathway_atlas/ENGINE_VS_LITERATURE.md` sec. 3 rows E6 and E11, and the binding red-team amendment,
sec. 5 item 2. Amended the same night by the plan premortem (coordinator, binding):
(1) the cells that diagnosed E6 (Cerny pot D, Givaudan, Hofmann) are development checks; held-out tests come from the
labelling cells of CAAS / BTBU;
(2) the change is scored at the incumbent constants first, and a compensation ledger replaces the blanket B60 guard;
(3) the refit runs only after the thiol leave-one-lab-out (EX-056) has landed.
The eventual target is ONE joint refit after E1 + E6 + E11. This wave is the E6 + E11 step of it, not a substitute.
The baselines in sec. 7 are the SHIPPED engine's, from read-only probes (sec. 10). No changed configuration has been
run.*

## 0. Summary for the morning

- **Structural finding (from the code):** the existing B10 map cannot make the route mix move with temperature.
  `FORMATION_ROUTE_OF` (`parameters_sulfur.py` ~l. 3141) puts `k_dpo_nf` and `k_dpo_ddp` both in `sugar_trunk`, and
  `k_nf_mft` and `k_ddp_mft` both in `thiol_assembly`. Each MFT route is one sugar step plus one thiol step, so two route
  barriers give both routes the same composite barrier, and the norfuraneol (NF) share stays T-invariant. The
  per-*step* barriers go through the B62 hook instead: `FurfuralBarriers.barriers_kj_mol` takes any fitted key, and
  `with_fitted_sulfur` reads it via `measured_ea_overrides()`. A generator sets this in memory, exactly as
  `generate_kinetic_core_b62_fit.py` does, so no src edit is needed for the fit. A SHIP would add a named switch.
- **Two of the red team's three numbers are design errors:**
  - **Hofmann 145 C.** On the actual bundles (NF + Cys 20 / 20 mM in 50 mL vs ribose + Cys 100 / 33 mM in 100 mL;
    ppb basis), the shipped NF / ribose MFT ratio is **11.5 vs 5.13 measured (+0.35 dex): it passes**. The red team ran
    both pots at 20 mM and compared with a per-pot microgram ratio.
  - **Givaudan.** The pot is ~2 M (37 % water), not 300 mM. The verdict does not change (0.0026 vs >= 48).
- **Held-out power is low:**
  - Every laboratory that measures the free-NF share diagnosed E6.
  - The CAAS / BTBU labelling cells test the intact-C5 vs fragment share. Both MFT routes are intact C5, so these cells
    are side-effect guards with ~0 power on E6's claim.
  - One held-out cell has leverage: BTBU's Wang 2012 xylose share of MFT (shipped 0.68 vs <= ~0.02; 68 % of the
    shipped MFT there runs through NF).
- **Expected verdict: INCONCLUSIVE or DO NOT SHIP (P 0.9).** SHIP-eligible P ~0.05. Compute: ~2.5 - 5.5 h wall on
  8 cores (worst ~11 h), launched after EX-056.

## 1. The change

**Base:** the shipped wave, B60: its 52-row objective, physics, numerics, declared sigmas, bands, optimiser and
two-start protocol (`generate_kinetic_core_b60_fit.py` via the B62 generator pattern). If B64 has shipped by launch,
its objective and optimum are the base, and sec. 7 is re-run first.

**Unfrozen (2 coordinates, added to B60's 23 free):**
- `k_dpo_nf`: log10 k at 145 C, B60 frozen +0.434, bounds (-10, 0.5) (`FITTED_SULFUR_BOUNDS_LOG10K`).
- **the DDP-route multiplier**: one log10 offset on `k_ddp_mft` and `k_ddp_mft_hs` together. B60 frozen -6.541 / -7.693.
  The pair is tied as in the joint fit's `DDP_COORDINATE`; bounds those of each key. `k_dpo_ddp` stays frozen (B60
  +0.042), so the route has one free number. This is the identifiable product, as the round-3 review found.

**Barriers (declared from measured classes; NOT fitted; all through the B62 hook):**

| step | B60 | this wave | class bound and source |
| --- | --- | --- | --- |
| `k_tdp_fur` | 64.08 | **147.0** | furfural formation, Schirle-Keller & Reineccius 1992 Table 2 (B62) |
| `k_fur_decay` | family 143.5 | **46.5** | furfural loss with Cys, Yaghmur 2005 (B62) |
| `k_dpo_nf` | 64.08 | **83.8** (sensitivity 93.0) | NF formation 81 - 93: Lau 2003 81.5 +- 3.2, Pandit 2006 93.0 +- 6.4 (`lau2003_extraction.md`, `pandit2006_extraction.md`); 83.8 is the inverse-variance mean |
| `k_nf_mft` | 64.08 | **profiled {64.08, 100, 133}** | thiol-assembly class [40, 133] (thiol prereg sec. 13); Chan & Reineccius 1994 sulfur volatiles 81 - 137. Profiled, not freed: EX-030 found it unidentified (cost flat within 4.8) |
| `k_dpo_ddp`, `k_ddp_mft(_hs)` | 64.08 | 64.08 (unchanged) | no measurement exists; stays on B60's lumped value (rule e: not a hiding place, it is not freed) |
| H2S supply `k_cys_h2s` | 133.0 | 133.0 (unchanged) | Cys thermolysis 123.3 - 135.5, Zheng & Ho 1994 Table V; the class bound if any H2S-release step is touched (none is) |

k(145 C) of each re-barriered step is held at the fit's log10 k, so a barrier change moves only the off-145 C rates.

**Composite barrier gap, NF route minus DDP route:** +19.7 / +55.7 / +88.7 kJ/mol at the three `k_nf_mft` points.
From 95 to 145 C that multiplies the NF / DDP odds by x2.2 / x8.8 / x32. That is the T-dependence the single 64.08
cannot produce.

**Members:**
- N1: B60 + B62 only, the B60 free set (B62's own candidate A, so B62's effect is separable);
- A64, A100, A133: the full change at each `k_nf_mft` point, start 0 (B60 optimum, the two new coordinates at their
  shipped values);
- A*-s1: start 1 (B8's perturbation protocol) of the lowest-cost converged A;
- A*-93: `k_dpo_nf` at 93.0 for the same member.

## 2. Stage 0: incumbent constants first (no refit)

Before any fit, run B60's frozen vector with the barrier change and without re-fitting: B62 on, `k_dpo_nf` 83.8,
`k_nf_mft` at each profile point (3 configurations, plus B62 alone). Score every test in sec. 4 and all 52 B60 rows.
This records what the barriers do to the shipped constants; the refit (stage 1) then shows what compensation buys.
A stage-1 row that gets worse is judged against stage 0 *and* B60.

## 3. Rows and measurement model

- **Rows:** B60's 52 (the B60 fit report's `row_ids`; `REMOVED_ROW_IDS`, `ROW_CORRECTIONS` and
  `SYSTEM_CORRECTIONS` as B60).
- **Measurement model:** B60's declared per-row sigmas, FIXED and identical across every member. This is EX-024's
  principle: one fixed model, so costs compare across members; no noise is re-estimated.
- **Converged-fit rules (EX-024):** eval budget 1200; ftol 1e-4; only status 2 / 3 counts. A budget-stopped member is
  not read, and is re-launched warm from its last point once. The verdict reads the lowest-cost converged A, and
  start 1 must agree within 2 cost units, or the verdict is INCONCLUSIVE.
- **Not in the objective:** every test in sec. 4.

## 4. Tests (frozen now)

### 4.1 Held-out (labs not cited for E6's diagnosis; binding where powered)

- **H1, CAAS: Liu 2023b, 168 C (CF-086).** Ribose 19.6 mg + Cys 198 mg in 0.5 M K3PO4 pH 5.88, 20 - 60 min; volume
  unstated, declared 10 mL, swept 5 / 20.
  - Test: C2 + C3 fragment share of MFT (`r_ha_mp_mft`) <= 0.12 under the corrected reading (~6 %, the impurity level;
    P 0.8). Under the printed reading (19 %) the pass band is [0.10, 0.38].
  - Shipped: ~1e-13, which passes the corrected reading and fails the printed one.
  - **Power against E6 / E11: ~0.** No member touches the fragment route (`r_ha_mp_mft`, thiol assembly at 64.08), so
    every rival predicts ~0. H1 is a guard: it fails only if a member opens a fragment route.
- **H2, BTBU: Wang 2012 pot C.** Xylose 26 mM (half [13C5]) + Cys 12.1 + Gly 12.1 + thiamine 1.19 mM, 25 mL,
  unbuffered, pH 4.8, 132 C, 90 min.
  - Test: xylose share of MFT (1 - `r_hmp_mft` share) **<= 0.04** (measured < ~0.02, ">99 % non-xylose", 0.3 dex
    tolerance). Pots A / B (GSH; Glu + Cys + Gly) are scored with the same rule when the engine answers them.
  - Shipped: **0.68 (fail, +1.2 dex)**; NF route 0.68 of MFT, DDP 1e-4.
  - The **only held-out test with leverage on the NF route.** It is confounded: the thiamine route's constant was
    fitted on Cerny's ternary share (B60 row `cerny_ternary_thiamine_share`), so a pass could come from the thiamine
    side. Read the member's thiamine flux beside it.
- **H3 / H4, BTBU / CAAS hexose labelling: Zhao 2019b (GSH + glucose, 140 C) and Du 2023 (Cys + glucose, 150 C).**
  Intact glucose C5, mixed channels < 1 - 6 %. **Not scored:** the shipped hexose MFT is 100 % C2 + C3, and no member
  can change that. Attributed to E1, which is the joint refit's job.

**Power statement.** The primary claim is that the free-NF share is low at <= 100 C and higher at 145 C. No held-out
laboratory measures that share:
- Nestle / Firmenich (Cerny), Givaudan and DFA diagnosed it;
- Reading's fed-NF pots (Whitfield 1999 / 2001) are B60 fit rows.

H2 is the one held-out row the claim moves. P(it discriminates the NF-route change from thiamine-side compensation)
~0.4 (mine). This wave cannot *confirm* E6 on held-out data. It can refute it (H2 or the ledger failing) or leave it
INCONCLUSIVE.

### 4.2 Development checks (the cells that diagnosed E6; reported, gating as stated in sec. 5)

- **D1, Cerny 2003 pot D (Nestle).** NF 145 + [13C5]ribose 145 + Cys 96.6 mM nominal, 0.5 M K phosphate pH 5, 95 C,
  240 min.
  - Test: fed-NF MFT share (the joint fit's tracer, `fed_norfuraneol_tracer`) in [0.03, 0.15] (measured 0.07).
  - Shipped 0.981.
- **D2, Givaudan (de Roos 2005; Laarhoven, unpublished).** Cys : carbonyl : lactate 1 : 1 : 1.67, 37 % water (~2 M
  declared; 0.3 M sensitivity), pH 5 clamped, 100 C, 90 min.
  - Test: ribose / NF MFT **>= 10** (measured 48).
  - Shipped 0.0026 (2 M) / 0.0066 (0.3 M).
- **D3, Hofmann 1998, 145 C / 20 min / pH 5 / 0.5 M phosphate, on the bundles through the panel path.**
  - Test: MFT(NF + Cys) / MFT(ribose + Cys) in ppb within 3x of **5.13** (that is, [1.71, 15.4]).
  - Shipped **11.5 (pass)**. Direct integration without the vessel oxidant gives 59.9: the panel's vessel-O2 model
    matters, so the test uses the panel path only.
- **D4, F23#13 (DFA): in-situ NF-route share of MFT in Hofmann's ribose + Cys pot, 145 C.**
  - Test: share in [0.30, 1.0].
  - Shipped 0.9998 (pass). A single lab's unlabelled in-house estimate (`thiols.md`), contradicted at the same pH by
    Cerny 2007's NF marker being absent at pH 4 - 5.5, 145 C (F23#10).
  - **D1 with D4 is the T-dependence signature** the single barrier cannot produce.
- **D5, Cerny pot C (ribose + Cys, 95 C).** MP3P / MFT <= 0.1 (MP3P not detected). Shipped 0.075.

### 4.3 Guard (compensation ledger; replaces "no B60 row worse by > 0.3 dex")

A B60 row may get worse than its B60 raw |log10 residual| by > 0.3 dex (in stage 0 or stage 1) only if it is listed
here, with its named open error, before the run, and only in the direction of its shipped residual.

| row(s) | named error | why the change exposes it |
| --- | --- | --- |
| `zhou_pH7_MFT`, `zhou_pH7_FFT`, `zhou_pH7_MFT_over_FFT`, `zhou_pH7_FUR_arp_alone` | **E5** (no buffer catalysis: an unbuffered ARP + Cys pot on constants measured in 0.5 M phosphate) and **E2** (no pH on the sulfur lane's trunk steps) | they sit at 120 C, off the 145 C reference, where the barrier change acts; shipped residuals already -1.4 to -1.9 |

Every other row is binding, explicitly including `kang_*` (100 / 120 C) and `zhai_*_fold_120_over_100`. Those rows
*are* the temperature dependence the change claims to fix, so a regression there counts against it. E1 and E7 have
no row in B60. A row not in the ledger that gets worse by > 0.3 dex means DO NOT SHIP.

## 5. Decision rule (frozen)

- **SHIP-eligible** (owner approval: every pentose MFT and the MFT / FFT ratios move):
  - a converged A member with H2 pass, H1 pass and the ledger guard holding;
  - D1 AND D2 pass, so the development checks do not contradict it;
  - D3 holds.
  - D4 is reported. A D4 failure blocks only if D1 passes by a share < 0.03 (the two-sided band).
- **INCONCLUSIVE:** D1 and D2 pass and the guard holds, but H2 fails, or H2 passes only through the thiamine route
  (thiamine MFT flux x > 3 vs B60 in that pot). Then E6 stands as a development-only result, and it is folded into the
  joint E1 + E6 + E11 refit with H2 and the CAAS / BTBU cells as that refit's held-out set.
- **DO NOT SHIP:** the guard fails at every A member, or no A member converges.
- **REFUSE (declared band):** D1 fails at every `k_nf_mft` point after the refit. The structure (two intact-C5 routes,
  class barriers) then cannot carry Cerny's 95 C route mix. The free-NF share of MFT is declared unidentified in
  95 - 145 C, the envelope carries MFT through both routes at [0.03, 1], and the designed experiment goes to the brief:
  Cerny pot D at 120 and 145 C, plus 2-mercapto-3-pentanone (sec. 3.1, E6 / E11).

## 6. Expectations (mine; informed by EX-024 / 030 / 035 / 036 / 041, where pot D failed 0.98 - 0.99 in every joint arm, 0.805 at Ea 133, and by EX-040's diagnostic pass of 0.070 at Ea 133 with DDP +3, which was not a fit)

| statement | P |
| --- | --- |
| Stage 0 passes D1 at any `k_nf_mft` point (DDP still ~4 decades weak) | 0.05 |
| Every A member converges within 1200 evaluations | 0.7 |
| The refit raises the DDP route by >= 2 decades at A133 | 0.6 |
| D1 passes at the best A | 0.25 |
| D2 passes at the best A | 0.15 |
| D3 still holds | 0.8 |
| D4 holds, given D1 passes | 0.5 |
| H2 passes | 0.15 |
| The ledger guard holds | 0.45 |
| Verdict SHIP-eligible / INCONCLUSIVE / DO NOT SHIP / REFUSE | 0.05 / 0.35 / 0.35 / 0.25 |

## 7. Shipped baselines (probe, sec. 10)

| test | shipped B60 | target |
| --- | --- | --- |
| H1 Liu 168 C fragment share (5 / 10 / 20 mL) | ~1e-13 (NF share 0.96 / 0.99 / 0.996) | <= 0.12 |
| H2 Wang 2012 pot C xylose share of MFT | 0.68 | <= 0.04 |
| D1 Cerny pot D fed-NF share | 0.981 (route flux share 0.9999) | [0.03, 0.15] |
| D2 Givaudan ribose / NF MFT | 0.0026 (2 M), 0.0066 (0.3 M) | >= 10 |
| D3 Hofmann NF + Cys / ribose + Cys MFT (panel) | 11.5 (direct, no oxidant: 59.9) | [1.71, 15.4] |
| D4 Hofmann ribose 145 C in-situ NF share | 0.9998 (100 C 0.9996, 120 C 0.9996) | [0.30, 1.0] |
| D5 Cerny pot C MP3P / MFT | 0.075 | <= 0.1 |
| B60 cost (52 rows, quick setting, re-evaluated now) | 27.132 (report 27.13) | ledger rule |

All route barriers today: `k_dpo_nf`, `k_dpo_ddp`, `k_nf_mft`, `k_ddp_mft(_hs)`, `k_pent_dpo`, `k_tdp_fur` 64.08;
`k_nf_decay` 143.5; `k_cys_h2s` 133.0.

## 8. Compute estimate (launch after EX-056 lands)

- **One B60-objective evaluation: 23 s serial** (measured now, twice, at the shipped vector). Earlier members took 246
  evaluations (B60) and 495 - 600 (B64), so assume 300 - 600 per member under the 1200 budget.
- **Per member:** 1.9 - 3.8 h serial; ~0.8 - 1.8 h at 4 workers (M2's 2.1 - 2.5x).
- **6 members** (N1, A64, A100, A133, A*-s1, A*-93), run two at a time x 4 workers (8 cores, the cap): 3 waves, so
  **~2.5 - 5.5 h wall**. Worst case, members budget-stopped at 1200: ~11 h.
- **Before launch** (work, not compute): ~2 - 3 h.
  - The generator, modelled on `generate_kinetic_core_b62_fit.py`, plus the route-test scorer ported from the probe.
  - V-checks: with the extra barriers set to 64.08 and the new offsets at shipped values, the cost is 27.132 to 1e-9;
    the scorer reproduces sec. 7 to 1e-6.
  - Stage 0 (4 configurations x ~30 s), ~5 min.
- If B64 has shipped: ~13 s / evaluation (inline coupling), so ~40 % less.

## 9. How it could mislead

- (a) **A free constant hiding the answer.** The DDP route has no measured barrier or rate. If D1 / D2 pass by raising
  it ~4 decades at 64.08, the route mix is right at 95 - 100 C for a reason the data cannot check. Rule (e) is only
  half met.
- (b) **Wrong temperature shape.** D4's 145 C support is one lab's unlabelled estimate, contradicted at pH 5 by Cerny
  2007. A member that "passes the T-dependence" may be fitting the wrong shape.
- (c) **Concentration.** Cerny's nominal concentrations (a ~0.5 g charge) and Givaudan's ~2 M at 37 % water are
  extrapolations, and the fed share is concentration-dependent.
- (d) **Path dependence.** D3 is path-dependent (panel vessel O2 vs direct integration: x5). Score it on one path only.
- (e) **H2's pass may be the thiamine route's** (sec. 4.1).
- (f) **Sloppy directions.** The two new coordinates trade against `k_nf_mft` and the FFT branch. Read the Laplace
  before claiming identification.
- (g) **This lane alone.** The joint fit (R3En) already frees both offsets on a different engine state (B58 / B59b /
  B59a / B61 on) and fails D1. A pass here on the B60 state need not survive the eventual joint E1 + E6 + E11 refit.

## 10. Probes (read-only, one core, nice 10, in memory; nothing written to the repo)

All in `/private/tmp/claude-501/-Users-pabloantoniomorenocasares-Developer-Maillard/128e3455-67ec-4d82-ba48-6cfcfb79a63e/scratchpad/prereg_probes/`:
- `probe_b.py`: D1, D2, D4, D5 and the barriers;
- `probe_b_panel.py`: D3 through `scoring.score_panel`;
- `probe_b_heldout.py`: H1 and H2;
- `probe_b60_eval.py`: the 23 s evaluation and the cost 27.132.

Each writes a `.json` / `.log` beside itself.
