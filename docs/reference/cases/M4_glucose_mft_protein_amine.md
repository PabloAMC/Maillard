# Case M4: glucose to MFT; protein and other amines as classes

**FROZEN 2026-11-06: items 1 - 5. Development grade on disk; decider W6 (decisive-measurements entry
`glucose_mft_carbon_origin`), sized: arm 2 + arm 3, n = 3, 96 heatings.**
- **Freeze commit:** **c7c2bf2537cf384d701267f67d913b3a9c0a533d** (c7c2bf25). It is the commit that froze items 1 - 5.
  A commit cannot contain its own sha, so this line was written by the commit that follows it. That commit changes
  nothing else.
- **Status of each step:**
  - the reduced-model power check has run (sec. 5.4 results, PR #101) and its demotions are written in;
  - the premortem was adopted (FM1 - FM9);
  - the coordinator's rulings are applied.
- **From here:** items 1 - 5 change only by amendment, declared before any new data are read. Items 6 - 9 are appended
  later, never edited into this part.
- **Open on the freeze:**
  - the power numbers must be re-checked under the revised judge once PR #99 is adopted (sec. 5.4, "a rule problem");
  - the W6 request is sized from them.
- **Item 5 tolerances follow sec. 7 item 5 as revised** (branch `judge-revision`, commit 455d8396; adoption pending its
  confirmation).
  - Before any scoring, the case's measurement model (sec. 5.1) and the reduced power check (sec. 5.4) are re-checked
    under the adopted text.
  - The revised judge also removes the "means differ by more than omega" adequacy escape that the power check flagged.
  - **No decisive data is read before then.**

No engine was run for this case, on any pot. The changelog is at the end.

**Coordinator rulings (2026-11-06), as corrected by the adopted premortem (FM1 - FM9).**
1. **The power check runs on a reduced forced-input mini-model** (sec. 5.4, ~1 day, M3-layer style), not on the full
   network. The full network (W3) is still needed for the M4 build itself. *(Supersedes "the power check waits for W3".)*
2. **The spent-data rule, as applied here (FM2): spent or seen data may FALSIFY a rival whose prediction there is
   constant-free (structural), and never CREDIT one.**
   - **D3** (the BTBU / CAAS label shares, F23#7 / #8; ENGINE_VS_LITERATURE l. 259 / 298) is what nominated E1, so it
     is spent.
   - **D1** was seen through R0's no-refit probe.
   - **D2's pH 3 / 7 arms** are MF-12 observations.
   - So none of the three can credit a rival. A structural prediction they contradict is falsified:
     - **C2 + C3 as the hexose MFT route is falsified**, both as R0's route (`r_glc_c2c3`) and as R3, by D3's
       constant-free label share;
     - **R4-as-written is falsified** by D5 (no pool).
   - **No rival is credited on disk.** Part G is **development grade**, and the decider is the designed experiment
     (sec. 4.4; decisive-measurements entry `glucose_mft_carbon_origin`).
   - D1 and D2 are scored and reported as development units.
   - *(This corrects the earlier basis, "D1 / D2 are the misfit that nominated the claim".)*
3. **The older hold-out declaration binds against fitting.** No rival's constants may be fitted on the glucose pH 3 / 7
   bundles; they are scored only (falsify-only under ruling 2).
4. **The claim is restated** as "a hexose C5 entry into 1-DPO; carbon origin (C-1 retained vs lost) undecided"
   (sec. 1). R1 (C-1 retained, through 1-DG) vs R5 (C-1 lost, through glucosone) is an explicit pair that **cannot be
   decided on disk**.
5. **Premortem additions:**
   - **FM3:** rival R4' (sec. 3).
   - **FM4 / FM5:** temperature-dependent mixtures, O2 on 1-DG, and the label / OPD / standard fixes (secs. 4.4, 5.1,
     5.4).
   - **FM6:** R5's band (sec. 5.2).
   - **FM7:** part A is an engine-structure requirement only (sec. 1).
   - **FM8:** **no M4 comparison runs before M2's verdict**.
   - **FM9:** the reservoir notes and one Cys-sugar Amadori class module (secs. 5.2, 6).

*Plan step M4 ([STATUS_AND_PLAN sec. 5, W2](../../guides/STATUS_AND_PLAN.md)). It covers:*
- *errors E1 (F23-a, F24-a) and E8 (F25-a, -b, F06-b) of
  [ENGINE_VS_LITERATURE](../pathway_atlas/ENGINE_VS_LITERATURE.md) sec. 3, as amended in sec. 5 (E1: "engine error
  SUPPORTED; band-now fix REFUTED");*
- *the plan item "amine as a variable" (F02-b, F05-a, F19-e);*
- *misfit MF-12 of the M0 inventory and its group G8 (sec. 2.11).*

*Items 1 - 5 of sec. 7. This text starts from the 2026-11-05 refresh (branch `m4-case-refresh`, 8660c005, not
merged), which is newer than the 2026-11-04 text on main. Items 6 - 9 are appended later and never edited into this
part. M4 is built on the one network (W3, [A2_MODULE_CONTRACT](../../guides/A2_MODULE_CONTRACT.md)); sec. 6 says what it
needs from it. Part E (E8) is kept for the record but is not built: the plan's M4 row defers it until M1's
cysteine-Amadori route is decided.*

**Rules this case applies (M1.1, 2026-11-05).**
- **Value basis per series.** Each decisive series states its basis (sec. 5.1): an absolute level, a within-pot ratio, a
  within-study contrast, a C/C0 course or a label share. H0's between-laboratory tau (tau_ABS 0.964 dex) enters only
  on absolute levels compared across laboratories or with an engine level. It never enters on a ratio, a contrast, a
  C/C0 course or a share.
- **omega = 0.419 dex** (= 0.965 ln) per printed value in every unit's primary covariance; 0.390 is a sensitivity.
- **Floors** come from each assay's reported LOD / LOQ or reporting limit, never from a data point. Where a paper
  prints none, the floor is a declared transfer with required agreement (sec. 5.3).
- **Densities** are exact, or Monte-Carlo with PSIS k-hat < 0.7 and ESS >= 100, on every form that can flip a sign.
  That includes censored (below-floor) forms and every sensitivity, not only the primary.
- **Adequacy.** A unit favours a rival only where that rival is adequate (posterior-predictive p > 0.05 under the
  unit's own measurement model) or the rivals' means differ by more than omega. Otherwise the unit reads "both
  inadequate". Holm across guard groups.
- **The pinned incumbent.** A favoured rival must also beat R0 pinned at its shipped (B60) constants, not only R0
  refitted.
- **Bands and starts** never come from data scored or screened in this case (secs. 4.1, 5.2). An unsourced number
  enters with required agreement: the verdict must hold at both ends of its declared range.
- **First-reader band sources.** A band whose only source is a first-reader dossier is carried as a required-agreement
  sensitivity until its second reader. The verdict must hold both at the dossier's band and at the band this case had
  before that dossier.
- **Screen first, demotions only.** Before any fit there is a mechanism-signature screen (sec. 5.5). After it a rival
  can be demoted, never promoted or redesigned, and the screen's pots are declared seen.
- **Power check** before the freeze (sec. 5.4).
- **Spent data** (sec. 7 item 2): data that an M0 check or any earlier diagnostic used to nominate a rival are
  development evidence for that rival and cannot be its decisive data (sec. 2.10).
- **pH scale.** Unless stated, every pH below is a **label pH**: the buffer or solution as prepared, at room
  temperature, before heating.
  - "Held" = a pH-stat, or a buffer with a reported drift.
  - "Cooled" = read after heating and cooling.
  - "In situ" = at temperature (an engine quantity unless a paper says otherwise).
  - The engine is compared with a measurement only on the measurement's scale (`observations.ph_cooled` for label or
    cooled readings).

**Split notation.** "v1" = [split.yml](../contrast_benchmark/split.yml); "v2" =
[split_v2_per_change.yml](../contrast_benchmark/split_v2_per_change.yml).
- v2 E1: DEV labs BTBU, CAAS, DFA Garching, Nestle / Cerny; **no TEST contrast** (UNTESTABLE_on_held_out).
- v2 E8: DEV labs DFA, Reading, Wageningen; **no TEST contrast**. Cerny 2007's F13#7, F25#6, F25#7 and F25#22 are
  excluded as calibration input.
- Atlas cells are 0-based indices of the F files in [pathway_atlas](../pathway_atlas/). Every cell below is v1
  DEVELOPMENT unless marked.

Two cautions:
- F23#7, F23#8, F24#4, F24#5 and F25#23 (Zhao 2019b, Du 2023) are DEV(E1) but **TEST for E6 and E11 in v2** (M2's
  errors). Their keyed tables (`heldout_du2023_tab2.yml`, `heldout_zhao2019b_tab1_2.yml`) are TEST-only
  ([contrast-benchmark README](../contrast_benchmark/README.md) sec. 11).
  - This case reads only their route pattern (intact C5 vs mixed isotopomers).
  - Anything they say about the norfuraneol share, or about the route mix with temperature, is M2's held-out evidence.
    It is read in item 7 only.
  - **No engine run on them is allowed before M2's item 7.**
- The Hofmann glucose pH 3 / 7 arms are DEV under both frozen splits. They are also, separately, declared Maillard-path
  hold-out bundles (`data/benchmarks/external_validation/maillard_path/mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH{3,7}.json`).
  The pH 5 arm is a trust-loop predictive row, not a B60 fit row (`kinetic_core_b60_fit_targets.json` holds only the
  ribose and xylose T1 rows). **Ruling 3:** the hold-out declaration binds against fitting. No rival's constants are
  fitted on the pH 3 / 7 bundles; they are scored only, falsify-only for part G (ruling 2).

## 0. What has been seen under the engine

**By the authors of this case.** Engine outputs read:
- ENGINE_VS_LITERATURE Appendix A, probes C and E: glucose + Cys MFT 1.5e-14 mM at pH 6, 100 % C2 + C3; ribose /
  glucose MFT 4.4e10 and FFT 0.47 on the Hofmann T1 design, pH 5.
- Appendix A, probes B and D: ribose + Cys FFT and MFT vs pH 3 - 7, 30 mM.
- The sec. 5 red-team probe of `HEXOSE_C5_ENTRY` (a scratchpad probe, not an artifact): ribose / glucose MFT 45 - 78
  without amine, 0.79 with 33 mM glycine, T1 design.
- `core_panel_scores.md` rows for the three Hofmann glucose bundles (all "refused", no number) and for
  `thiamine_cys_glucose_120C_Bolton1994` (predicted 232.9 vs 11.7 ppb).
- **Read while refreshing (2026-11-05):** the per-data-set summary of the no-refit thiol sub-model probe
  (`kinetic_core_thiol_submodel_probe.md`, shipped constants = R0). Its rows:
  - Hofmann Table 5 NF: mean residual -12.3 dex over 2 rows (worst -24.8);
  - Table 2 MFT: -6.0 over 4 rows (worst -9.7);
  - Table 1 MFT: -5.4 over 4.

  So D1's and D2's pots have been run under R0 and that summary is seen. No per-row number was read, and no rival
  other than R0 has run on them.
- **Read while preparing the freeze (2026-11-06):** the M0 misfit inventory (draft PR #82, branch `misfit-inventory`).
  - Its engine residuals on G4 and G8 members (sec. 2.11): MF-12 (the probe numbers above), MF-15 (Bolton 1994
    +1.30 dex), MF-23 (Steinhagen 2021: glucosone -1.80, glyoxal -1.54, HMF -0.97, 3,4-DDG -0.82, 3-DG +0.10,
    MGO 0.00 dex), MF-24, MF-29 / MF-30 (M3's Strecker residuals) and MF-34 (pyrazines -2.9 dex).
  - None is on a pot this case uses as decisive.

No engine output exists on any v2 E1 or E8 TEST contrast (both sets are empty). Measured values have been read only for
the v1 TEST cells F06#16 (Lee 2022) and F14#5 (Hollnagel 2000), and for the v2 E5 TEST cell F06#11 (Thornalley 1999),
in the atlas files. Both sweeps (sec. 2.9) were screened on titles and abstracts only. The fourteen M4 sweep dossiers
(PR #76) were read in full by their first readers and by the authors of this case. They hold measured values only, no
engine output.

**Pots by engine status.**
- "Seen" = an engine number on that pot exists in a record (fit row, probe or screen), so the pot cannot be decisive
  here.
- "Run under R0 only" = the incumbent's no-refit probe has a number there and no other structure has run. The
  2026-11-05 refresh kept such pots decisive between the new rivals. **Ruling 2 overrides that:** seen data
  (D1 through R0's probe; D2's pH 3 / 7 arms as MF-12 observations) may falsify a structural prediction but never
  credit a rival. They are scored as development units.

| pot (record) | status | where |
|---|---|---|
| Hofmann 1998 Table 1 (= Table 2's pH 5 arms) | **seen** | ribose, xylose: B60 fit rows; fructose, glucose: `thiol_submodel.yml` rows and B8's level rows; T1 design in probe E and the red-team probe |
| Hofmann 1998 Table 2, ribose pH 3 / 7 | **seen** | `thiol_submodel.yml` rows; probes B and D in shape (30 mM design) |
| Hofmann 1998 Table 2, **glucose pH 3 / 7** | **run under R0 only; spent** (ruling 2: the E1 misfit; never fitted, ruling 3) | mp hold-out bundles; `thiol_submodel.yml` rows `hofmann1998_cys_glucose_ph{3,7}` |
| Hofmann 1998 Table 5, ribose NF and furfural | **seen** | B60 fit rows |
| Hofmann 1998 Table 5, **glucose NF / furfural** | **run under R0 only; seen, falsify-only** (ruling 2); rhamnose not run | `hofmann1998_cys_glucose_ph5_fa_nf` |
| Hofmann 1998 Tables 3, 4, 8, 10 (fed ribose / NF / furfural + H2S, HA + MP) | **seen** | B60 fit rows (`hofmann1998_c2c3_recombination_*` = Table 8) and `thiol_submodel.yml` rows |
| Hofmann 1998 Tables 6, 7 (HA + MGO -> NF / furfural; MGO + H2S -> MP) | **not run** | -- |
| Hofmann 1998 Table 9 (1-mercapto-2-propanone, glucose + Cys, pH 3 / 5 / 7) | **seen** | `thiol_submodel.yml` rows; thiol sub-model probes |
| Zhang 2026 Tables S6 - S8 (GSH + glucose MFT, 130 - 150 C) | **seen** | `thiol_submodel.yml` rows; thiol sub-model probe |
| de Roos 2005 Table I, ribose and fed-NF arms; Fig. 4 | **seen** | red team E6 ("Givaudan ribose / NF"); B59b constant (Fig. 4) |
| de Roos 2005 Table I **glucose arm**; Table II | **not run** | -- |
| Du 2023, Zhao 2019b | **not run; reserved** (TEST for E6, E11) | -- |
| Meynier 1995 | **seen** (Cys pots in shape via probe B; glycine pots in M1b's case) | -- |
| Bolton 1993 / 1994 | **seen** | core panel |
| Li 2011 (glucose + Cys, free Cys loss, 100 - 130 C) | **seen** | B59 hold-out diagnosis (withdrawn all-free set); B59b prereg |
| Deng 2021; Davidek 2006a; Martins 2003b, 2005; Brands 2001, 2002c | **seen** | pH-acid sub-model probes; trunk fit rows |
| Xia 2022 | **seen** | B21 |
| Steinhagen 2021 (glucose alone, 121 C) | **seen** | core panel (MF-23) |
| Tressl 1994, 1999; Glomb 2010; Voigt 2009, 2010; Smuda 2010; Kanzler 2021, 2026; Pfeifer 2013; Du 2024b; Luo 2022, 2023; Cannon 2023; Ames 2001; Beck 1988; Nedvidek 1992 | **not run** (no engine species for most of their observables) | -- |

## 1. The claim

Three parts, built and judged separately on the one network.

| part | step(s) | today (B60, read from the code) | M4 |
|---|---|---|---|
| G. hexose entry to MFT | **restated 2026-11-06 (ruling 4): a hexose C5 entry into 1-DPO; carbon origin (C-1 retained vs lost) undecided.** New steps: Glc + Cys -> cysteine-hexose Amadori compound (D-fructose-L-cysteine), then into 1-DPO by either of two forms. **Form R1** (C-1 retained): -> 1-DG + Cys released -> `r_odg_dpo` (1-DG -> 1-DPO + C1, [sulfur.py](../../../src/kinetic_core/sulfur.py) ~927, "the C-6 loss"). **Form R5** (C-1 lost): -> glucosone + Cys (B21's `r_ama_g`) -> hydrolytic beta-cleavage (formic acid = C-1) -> 3-DPO / 1-DPO. Both then feed the shared pentose MFT steps | see below | a measured-class Cys-hexose Amadori formation and decay; the 1-DPO entry on with declared barriers, in both forms |
| E. the furfural supply's pH law (E8) | `r_pent_tdp`, `r_arp_tdp`, `r_glc_fur` (`acid` = [H+], a decade per unit; sulfur.py ~1160 - 1175) | as stated | an amine-catalysed 1,2-enolisation term, pH-flat or weakly rising, beside a weaker acid limb. **Deferred** (plan M4 row) |
| A. amine classes | the trunk's `r_schiff` (Glc + Gly), the sulfur lane's Cys-catalysed `r_pent_dpo` / `r_pent_tdp`, B20's `r_glc_lysp` (protein lysine, trunk only), the acrylamide lane's one-step `a_gln_glc` / `a_lys_glc` / `a_ala_glc` (M0a P1b) | see below | classes: free alpha-amino acid, cysteine, peptide, protein-bound lysine; each with its own condensation constant **and its own Amadori chain**; their competition for each sugar declared; **amine catalysis of dicarbonyl loss declared separately from amine consumption** (Kanzler 2026, sec. 2.2); **one pH module across lanes** (M0a P1 - P7) |

**Today, part G.** Glucose reaches MFT only through `r_glc_c2c3` (log10 k -8.59 at 145 C; Glc -> HA + MGO), then
`r_mgo_mp` (MGO + H2S -> 1-mercapto-2-propanone, MP) and `r_ha_mp_mft`. FFT forms through `r_glc_fur` (`acid`
factor). `HEXOSE_C5_ENTRY.k_per_min_100c = 0` ([parameters_sulfur.py](../../../src/kinetic_core/parameters_sulfur.py)
~623 - 641). Cysteine is not an amine on the trunk, so a Glc + Cys pot has no amine route to 1-DG.
`engine.unidentified_routes` refuses hexose-only thiol requests (engine.py ~1024 - 1043).

**Today, part A.** There is one amine: MET and PRO are folded into `Gly` (engine.py ~2560 - 2570). `LYSP` exists only
on the trunk lane. There is no peptide.

**Part A absorbs case M0a (parked 2026-11-06; draft PR #91, `docs/reference/cases/M0a_ph_law_every_lane.md` sec. 1.1,
CF-133).** M0a read the code and found two gaps that part A must close.
- **The trunk's pH law runs only in the trunk lane.** It covers the three Amadori-decay steps `k_ama_tdg` / `k_ama_odg` /
  `k_ama_mgo` (0.69 dec / pH) and the 3-DG -> formic exit `k_tdg_fa` (0.77 dec / pH), and is applied only
  `if lane == TRUNK`. Bypass paths, read from origin/main e560338e:
  - **P1:** the acrylamide lane never calls `trunk_conditions.apply`.
  - **P2:** in the sulfur lane `ph.trunk_terms` exists and is off.
  - **P3:** a pentose with Asn / Gln / Lys / Ala is refused as a lane conflict.
  - **P4:** lipid + acrylamide or sulfur inherits P1 / P2.
  - **P5:** the envelope draws the Amadori band for the trunk only.
  - **P6:** the JAX trunk module defaults to the acrylamide lane, with pH mode "none".
  - **P7:** `ph_final` is ignored outside the sulfur lane.
- **P1b: the competitor amines have no Amadori chain.** In the acrylamide lane Gln, Lys and Ala react only through a
  single amine + glucose -> melanoidin step (`a_gln_glc`, `a_lys_glc`, `a_ala_glc`; acrylamide.py l. 241 - 261). So no
  pH law on the Amadori steps can reach them; copying the law across lanes does not fix it.
- **Part A's requirement, therefore:**
  - amine classes as species **with their Amadori chains** (formation, decay to 1-DG / 3-DG / MGO, amine return) and
    their competition declared;
  - **one pH module across lanes**, so the Amadori-decay and 3-DG-exit laws act wherever those steps run, at the
    declared scale (initial, cooled or in situ);
  - P1 - P7 closed, or refused by name (sec. 6).

  M0a's rivals and decisive data (Schibilsky 2019's glucose + Ala DMHF pH response, spent there) are not re-argued
  here. Part A inherits M0a's spent list.
- **The pH scale, resolved (FM7).** The one pH module reads the **in-situ pH** wherever the lane carries a pH state (the
  sulfur lane, Amendment 1 of sub-model 1). Elsewhere it reads the **declared initial (label) pH**, held over the
  programme, as the trunk lane does today. P7 (`ph_final`) is honoured only where a programme declares it. Every
  comparison is made on the measurement's own scale (label, cooled or in situ; `observations.ph_cooled`). A lane that
  switches from label to in situ is a part-A case of its own, not a silent change.
- **M0a's sweep candidates, carried into part A with M0a's fixed roles** (unseen; none is read here):
  - Wang Y. 2019, 10.1016/j.foodres.2019.05.001 (HMF at >= 2 pH, Asn + hexose);
  - Wang P. 2022, 10.1002/jsfa.12149;
  - Ajandouz & Puigserver 1999, 10.1021/jf980928z;
  - Ajandouz 2001, 10.1111/j.1365-2621.2001.tb08213.x;
  - Chen 2011, 10.1016/j.foodres.2011.06.002;
  - Hollnagel & Kroh 1998, 10.1007/s002170050294.

  Each enters only by amendment, declared before its data are read.
- **Hollnagel & Kroh 1998, reconciled.** M0a lists it as a part-A evidence candidate: dicarbonyls at >= 2 pH with a
  non-glycine amine, from the Schibilsky laboratory, so it is not independent of the spent Schibilsky pots. This case's
  sweep 2 excluded it for **part G** under criteria (f) - (h): no thiol, no O2 contrast at >= 100 C, no fed C6. Both
  decisions stand. It is out of part G and a part-A evidence candidate under M0a's role.

**Why the claim is stated as a C5 entry with its carbon origin open (2026-11-06).** The 2026-11-05 claim named one
C5 step, 1-DG -> 1-DPO + C1, which keeps glucose C-1 as MFT's methyl. Two laboratories' positional labels speak
against that step, weakly:
- **Glomb 2010** (Halle; 50 C, pH 7.4, Lys, air vs de-aerated) finds the hexose C5 dicarbonyls, 1-DPO included,
  "strictly" glucose C2 - C6. They come from glucosone by hydrolytic beta-cleavage and rise with O2; none is assigned
  to 1-DG.
- **Tressl 1994** (TU Berlin; 160 C, pH 5, Cys) reports no C5 product from 1-DG; every 1-DG product keeps six carbons.

Both are absences without a floor, neither measures MFT, and neither is at Hofmann's 145 C, pH 3 - 7. So they do not
refute the C-6-loss form. But they name a second form with the same entry species and the opposite label (R5), and no
measurement on disk tells the two apart (sec. 4.3).

The claim is therefore what both forms share: the Cys-hexose Amadori pool feeds 1-DPO, and 1-DPO feeds the pentose
MFT steps. Which carbon is lost is left to the designed experiment (sec. 4.4, arms 1 and 2). R1 vs R5 is reported as
an explicit "cannot decide on disk" pair.

**What can ship from disk.** Nothing in part G (ruling 2): no rival is credited on disk, so R1 / R5 cannot ship over
R0. The on-disk comparison can only falsify structural predictions:
- **C2 + C3 as the hexose MFT route is falsified** by D3's constant-free label share, both as R0's route
  (`r_glc_c2c3`) and as R3;
- **R4-as-written is falsified** by D5.

It reports D1 / D2 as development units.

**Part A is an engine-structure requirement only (FM7).** It is built in W3 at parity with today's engine and off by
default. Switching it on (any amine-class constant moving a shipped output) is a separate part-A case with its own
items 1 - 9. Nothing in part A is decided in this case.

**Part E depends on M1.** Enabling `CYS_AMADORI_ROUTE` (M1's combination C1) retires `r_pent_dpo` / `r_pent_tdp`
(`inactive_reactions`, sulfur.py ~975 - 978). It also carries the pentose supply on de Roos 2005 Fig. 4's measured,
U-shaped Cys-Amadori decay. If C1 ships, E8's pentose part is decided by M1 / M1b, and part E keeps only `r_arp_tdp`
(fed ARP pots) and `r_glc_fur`.

**MGO on the one network.** Today the sulfur lane has no amine-free MGO sink.
- Its MGO comes from `r_glc_c2c3`, `r_dpo_c2c3`, `r_nf_decay_oh_mgo` and `r_tdp_mgo`. It leaves only by `r_mgo_mp`,
  `r_cys_actz` and the switched-off `r_mgo_cys_*`.
- On the one network MGO is one species. The trunk's `r_mgo_mel` (MGO -> melanoidin, first order, amine-free,
  **fitted**: log10 k(100 C) -1.643 +- 0.114, Ea **20.0 kJ/mol at its search floor**, `FITTED_EA_BOUNDS` (20, 260))
  then acts on every MGO the sulfur module makes.

The measurements, all from one laboratory lineage (TU Berlin food chemistry; CF-122, updated 2026-11-05):
- **Within one paper.** Kanzler 2026 measures an apparent Ea of **95 - 121 kJ/mol** for MGO + Ala / Phe / Pro, 80 ->
  100 C, pH 5 (initial, unbuffered).
- **MGO alone, across papers.** Against Bork 2022 (130 C), the estimate is ~96 - 114 kJ/mol with Kanzler 2026's 100 C
  point and **~47 - 66 kJ/mol with Kanzler 2021's**. Kanzler 2021 is the same laboratory and nominally the same pot,
  with a k2 3.2x higher (0.51 dex).
- **Amine-free at 130 C.** Pfeifer 2013 (same lineage) gives first order, k1 ~0.0053 /min, x7 below the shipped
  `k_mgo_mel` at 130 C. The shipped level is an amine-pot level.

From the same 100 C value, an Ea of 47 - 121 gives a loss x2.5 - 33 faster at 145 C than Ea 20 does [dry]. Every
rival that carries MGO inherits a sink whose temperature law the measurements contradict: R3 most of all, and R1, R2
and R5 through NF's MGO-releasing decay. Sec. 5.2 bands it; R0 is also scored pinned at Ea 20.

**Earlier waves on these structures** (cited from the record, none re-run):

| structure | wave | outcome | kind |
|---|---|---|---|
| hexose MFT through `k_glc_ha` (C2 + C3), carried by the eight Hofmann T1 level rows | B8 | the route existed only through the rows scored on it | whole-model fit |
| the same with the T1 level rows dropped | B9 (`tasks/data_restructure_plan.md`, 2026-09-03) | glucose and fructose MFT predict zero; `k_glc_ha` at its band floor; B9's Laplace covariance singular in that coordinate | whole-model fit |
| the hexose route reopened | B47, B51.1 | DO NOT SHIP: `k_glc_ha` +3.6 decades bought back one Whitfield row; directional claim SUG-09 flips | whole-model fit |
| hexose -> C5 entry (`HEXOSE_C5_ENTRY`) | sub-model 4 (backlog, unchecked) | never pre-registered or run | -- |
| the entry switched on alone (the "band-now" fix) | red team, ENGINE_VS_LITERATURE sec. 5; M0 CK-41 | ribose / glucose MFT 45 - 78 with no amine (measured 10.4), 0.79 with 33 mM glycine | read-only probe |
| protein-bound lysine glycation | B20 (SHIPS, trunk only) | Nguyen 2016's constants; answers only with `protein_g_per_l` | fit on printed rates |
| pentose + glycine entry | EX-012 / EX-013 | admitted on audit; norfuraneol sanity bound FAILED (25 - 36 mol %) | fit |
| the furfural supply's pH law | Appendix A probe B | FFT x44 from pH 4.5 to 6 and ~x400 to 6.5 against x6.1 measured (Meynier, held) | probe |
| the trunk's MGO sink | B1 (`k_mgo_mel`) | Ea at its 20 kJ/mol floor; measured 47 - 121 (CF-122, MF-24) | fit; registered conflict |
| the aqueous glucosone supply (`r_ama_g`, `r_g_go`) | B21 | ships; fitted on Hamzalioglu 2026's LacLys -> glucosone constants, licensed at pH 6.7 only, no O2 law (`ATMOSPHERE_FACTORS` = {air: 1.0}) | fit on printed rates |

**Errors stack (CF-107).** The B8 / B9 / B47 / B51.1 outcomes are whole-model fits and weak evidence about the hexose
route. Candidate compensating errors:
- E6 (MFT through free norfuraneol, 99.9 %) and E11 (one 64.08 kJ/mol formation barrier);
- the frozen pot-independent `k_mft_decay`, and the Strecker sulfur sent to an inert pool (SULFUR_FLUX_MAP sec. 6);
- MFT destroyed in excess H2S at B9's optimum (backlog, "next sulfur question");
- E2 (the enolisation split frozen in pH), and the trunk's MGO level (B21), with `k_mgo_mel`'s Ea at its floor.

The probe's 45 - 78 is an endpoint ratio that inherits E6 and E11, so it does not reject the entry either. Only the
labelling (sec. 2.1) is a measured share that makes a route inert: it rejects C2 + C3 as the hexose route.

## 2. All the evidence, by kind

### 2.1 Isotope labelling

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Du 2023 (F23#8, F24#5; ledger `du2023_table2_camola_isotopomers`) | Cys 1 mmol + 1:1 [13C6]/[12C6]-glucose, 0.2 M phosphate pH 6.5, 150 C, 120 min: MFT M+0 : M+5 = 52 : 48 (57 : 43); no mixed channel printed | hexose MFT is an intact glucose C5; C2 + C3 absent on the label basis | DEV(E1); TEST(E6, E11) |
| Zhao 2019b (F23#7, F24#4, F25#23; `zhao2019b_table2_camola_glucose_labelling`) | GSH + glucose 60 / 60 mM, phosphate pH 6.5 (molarity not stated), 140 C, 5 h: MFT ~100 % glucose C5; FFT ~94 % intact C5 | the same in a peptide-amine pot | DEV(E1); TEST(E6, E11) |
| **Tressl 1994** (TU Berlin, Tressl group; `tressl1994_extraction.md`) | **positional** [1-13C]- and [6-13C]-glucose + Cys, equimolar, pH 5, 160 C, 1.5 h. Table IV: the furfuryl skeleton (furfuryl alcohol, 2-hydroxymethylthiophene) and the thiolanones are unlabelled from [1-13C] and 100 % labelled from [6-13C], so they are glucose C2 - C6 (C-1 lost, assigned to 3-DG). The C6 furfurylthiols (19, 20) carry C-1. Every 1-DG product keeps all six carbons (Table V). **MFT, NF and glucose FFT not reported** | which carbon a hexose C5 loses. No 1-DG-derived C5 is reported (absence, no LOD) | DEV; a fourth laboratory |
| **Tressl & Rewicki 1999** (same laboratory; review; `tressl1999_extraction.md`; new here) | Fig. 14, no conditions and no citation: FFT and furfural from [1-13C]-glucose + Cys each ~50 % singly labelled. The authors read this as C-1 / C-6 equilibration on the 3-DG side, then C1 - C5 cleavage. One sentence says MFT's glucose pathway "was established", with no datum | **conflicts with Tressl 1994 on the furfuryl skeleton** (unlabelled there). Unresolved: neither share is used for the furfuryl skeleton. If equilibration is real, a positional label on a 3-DG-side C5 reads ~50 % whatever the cleavage (sec. 4.4, arm 1) | DEV; review-scheme value |
| **Glomb 2010** (Halle, Glomb group; `glomb2010_extraction.md`; new here) | glucose (unlabelled, [1-13C], [6-13C]) + Lys 42 mM each, 0.1 M phosphate pH 7.4, 50 C, 7 d, aerated vs de-aerated (Ar + DTPA), post-incubation OPD. The C5 dicarbonyls (1-DPO 0.03 / 0.01, 3-DPO 0.21 / 0.10 mmol/mol glucose, aerated / de-aerated) are "strictly" glucose C2 - C6, from glucosone by hydrolytic beta-cleavage (formic acid = C-1). 1-DPO rises x3 with air while 1-DG falls x4.4. No C5 product from 1-DG. Shares are asserted, not tabulated; no LOD | **a C-1-lost, methyl-bearing 1-DPO (methyl = glucose C-2) through an oxidative route: rival R5.** Against `r_odg_dpo`'s C-6 loss at 50 C (absence, no floor) | DEV; a fifth laboratory |
| Schieberle 2005 (DFA; review; `schieberle2005b_extraction.md`; new here) | glucose + Pro CAMOLA: HDMF intact under roasting (160 C, 10 % water); "about 60 %" from C3 + C3 in aqueous buffer, 145 C, 20 min (author's estimate). Its Table 3 (fed C3 + C3 -> HDMF, pH 3 / 5 / 7) repeats two cells of Hofmann 1998 exactly (Table 6's 364.5 ug, Table 9's 59.6 ug): **not used**; only the primary (Hofmann 1998 Table 6) enters | fragment recombination is real in water at Hofmann's 145 C for a C6 furanone with Pro. It bears on R3's premise by analogy only (no thiol) | DEV |
| Liu 2023b (F23#6, F24#6; CF-086) | ribose 168 C: fragments ~6 % after the [M-H]+ correction | fragments stay minor at high T (pentose) | DEV; TEST(E6, E11) |
| Cerny 2003 / 2004, 2007b (F23#0 - #3, F24#0 - #3) | pentose MFT intact C5; free NF 7 % (pot D) | the downstream steps M4 reuses | DEV |
| Bolton 1993 (F27#6; `bolton1993_tab2_mft_34s_split`) | 34S-Cys + thiamine in a glucose bouillon, a_w 0.83, 120 C: Cys supplies ~8 % of MFT sulfur; no MFT without thiamine (< ~2.3 % of the thiamine pots) | an upper bound on Cys + glucose MFT at 120 C in a matrix | DEV |
| Davidek 2005, 2006a (F07#5 - #7) | fed DFG and 1-DG: acetic acid ~70 % from glucose C1 / C2 | 1-DG's own exit (beta-cleavage), the competitor of any 1-DG -> C5 step | DEV |
| Voigt 2010 (Halle, Glomb group) | [13C]glucose (C-1, 2, 3, 5, 6) + Ala, 50 C, 7 d, air: lactic acid 48 % C1 - C3, 36 % C3 - C5, 8 + 8 % C4 - C6 | 1-DG's fragmentation pattern at low T | DEV |
| none | **no positional label on MFT or on NF from a hexose**, and no time-resolved hexose label | R1 / R2 (C-6 lost, methyl from C-1) against R5 (C-1 lost, methyl from C-2) | gap; the designed experiment (sec. 4.4, arm 1) |

### 2.2 Fed intermediates

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Hofmann 1998 Table 8 (F24#7, F24#8) | fed HA + **1-mercapto-2-propanone** (MP; not C2 + C3 + H2S, corrected 2026-11-05), 145 C: MFT 0.01 / 0.23 / 0.27 mol % at pH 3 / 5 / 7 (rising with pH) | the C2 + C3 route's last step and its pH shape (R3) | DEV; B60 fit rows (`hofmann1998_c2c3_recombination_*`) |
| Hofmann 1998 Table 6 | fed HA + MGO (1 mmol each, no sulfur), 145 C: **NF 0.13 / 0.78 / 3.18 mol %** and furfural 0.11 / 0.40 / 1.60 mol % at pH 3 / 5 / 7 | **C2 + C3 make NF too**, rising with pH: R3 is not "no NF" | DEV; not run |
| Hofmann 1998 Table 7 | MGO + H2S, pH 5: MP 1.8 / 4.0 mol % at 1 / 2 mmol H2S (super-linear) | the C3 sulfur step | DEV; not run |
| Hofmann 1998 Tables 3 / 4 (F13#2, F23#12, F25#4) | fed NF / furfural + H2S, pH 5 | downstream steps shared with ribose | DEV; fit rows |
| Deng 2021 (F03#11 - #12, F06#6, F07#11) | fed Met-glucose ARP, 100 C, unbuffered (pH measured only in the buffer arm): k 0.0048 (label pH 5.5), 0.012 /min (label pH 7.5); Met returned ~mole for mole | a hexose amino-acid ARP decay and its pH slope (the class analogue for the Cys-hexose ARP) | DEV |
| **Luo 2023** (Ho network; `luo2023_extraction.md`; new here) | purified Gly-, diGly-, triGly-glucose ARPs, ~3 mM, water, initial pH 3 - 9, 100 / 130 / 160 C, 30 - 180 min: Gly-ARP decay 1.5 - 8.5e-3 /min at 100 C; apparent Ea 34 - 56 kJ/mol; fast-then-slow, not first order; amine returned 0.54 - 1.0 | the hexose ARP decay band's Ea and Cys-return ends (sec. 5.2); peptide ARPs less stable at >= 130 C (part A) | DEV |
| Martins 2003b (F03#3 - #5, F06#5, F07#10); Davidek 2002, 2005 (F03#0 - #2, F02#11) | fed DFG: amine returned 71 - 96 %; split moves with pH | the hexose ARP steps the claim reuses | DEV; Martins a calibration input |
| Davidek 2006a (F07#5, #9) | fed amine-free 1-DG, 0.2 M phosphate: acetic 0.20 / 0.56 / 0.80 mol/mol at pH 5 / 6 / >= 6.5; 1-DG gone within 5 min at pH >= 6 | 1-DG's lifetime and its main exit | DEV |
| Voigt 2009, Smuda 2010 (Halle) | fed 1-DG + Ala or Nalpha-BOC-Lys, 42 mM each, 0.1 M phosphate pH 7.4, 37 / 50 C, air vs deaerated: exits are acetic, glyceric, lactic and formic acids (310 : 88 : 58 : 27 mmol/mol, 7 d, BOC-Lys), C4 dicarbonyls and their amides; negligible MGO; **no C5 product reported**; no amine-free arm | 1-DG's measured exits at low T; a ceiling on `r_odg_dpo` only through the unaccounted share (sec. 5.2) | DEV |
| Glomb 2010 (as sec. 2.1) | fed dicarbonyls + Lys, 50 C, pH 7.4: half-lives 1-DG 0.5 h, glucosone 8 h, 3-DG 40 h; fed glucosone gives C5 dicarbonyls "in the same ratio as in glucose reactions" (no numbers) | a ceiling on R5's glucosone -> C5 step and on `r_odg_dpo` (sec. 5.2) | DEV |
| **Kanzler 2026** (TU Berlin, Kanzler group; `kanzler2026_extraction.md`) | fed MGO 100 mM alone or 1:1 with ten amino acids (no Gly, no Cys), water, initial pH 5.0 (unbuffered, room T), 100 C, 0 - 120 min, and 60 - 140 C at 30 min. MGO alone: 65 / 31 % left at 30 / 120 min (second order). Amino acids speed MGO loss x2.6 - 9.7 at 30 min **while barely consumed** (Lys 98 % left, Ala ~100 %). Ala order 0.75 - 0.84; apparent Ea 95 - 121 kJ/mol | (i) the temperature law of the MGO sink every rival inherits; (ii) part A: amines catalyse dicarbonyl loss without being consumed; (iii) R3's MGO standing level at 145 C | DEV (no split cell); not representable as a pot (no Gly / Cys); a band source (sec. 5.2) |
| **Kanzler 2021** (same laboratory; `kanzler2021_extraction.md`; new here) | MGO 100 mM alone or + Ala / Lys 100 mM, pH 5.0 (unbuffered), 2.5 mL ampoules, 100 C, 0 - 300 min, triplicates: MGO alone 36.5 / 21.9 / 11.6 / 9.0 / 7.4 % left at 30 / 60 / 120 / 180 / 300 min, second order, k2 ~5.9e-4 /mM/min. That is **x3.2 Kanzler 2026's k2 in a nominally identical pot**. With Ala / Lys: 12.9 / 10.2 % left at 30 min | the MGO-alone Ea basis moves to 47 - 66 kJ/mol against Bork 2022; the within-laboratory between-paper spread of k2 is 0.51 dex; the order is identified (second) within one series | DEV; band source |
| **Pfeifer 2013** (TU Berlin, Kroh group, the same lineage: shared Haase, address and OPD method; `pfeifer2013_extraction.md`; new here) | MGO alone ~10 mM, initial pH 5, 130 C, aerated: first order k1 ~0.0053 /min (0.0042 - 0.0065), no replicates; pH 8: 13 % left at 30 min | the amine-free MGO loss at 130 C, x7 below the shipped `k_mgo_mel` there. **Not a second laboratory** | DEV; band check |

### 2.3 Trapping

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Beck 1988, Nedvidek 1992 (F06#0 - #2, F07#0 - #3) | in-situ OPD on hexose ARPs / glucose + beta-Ala: 1-DG : 3-DG rises x12.5 from pH 4.5 to 7; a secondary-amine ARP goes to 1-DG even at pH 4.5 | the hexose 2,3-enolisation's pH law; amine structure moves the split | DEV (one school) |
| Nedvidek 1992 (F23#16) | the 1,4-dideoxypentosone quinoxaline needs alpha-amino acids | a 1,4-dideoxy C5 made by amino-acid Strecker on a deoxyosone (rival R2) | DEV |
| Gobert 2009 (F06#3, F07#4) | 50 C: air lowers 1-DG x3.4 | O2 drains the 1-DG / ARP pool (low T only); consistent with R5's oxidative exit | DEV |
| Glomb 2010 Table I (as sec. 2.1) | twelve glucose dicarbonyls, aerated vs de-aerated, 50 C: glucosone 0.30 / 0.16, 1-DPO 0.03 / 0.01, 1-DG 0.05 / 0.22 mmol/mol | R5's O2 dependence at low T | DEV |
| Hollnagel 2000 (F14#5) | dry glucose + OPD: deoxyhexosuloses < 1 vs 11 mol % with glycine | amine catalysis of enolisation | **v1 TEST**: listed, not used |
| Thornalley 1999 (F06#11) | BOC-Lys raises 3-DG formation x6 at 37 C | amine catalysis | **v2 TEST(E5)**: listed, not used |

### 2.4 Measured intermediates and lifetimes

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Hofmann 1998 Table 5 (`hofmann1998_t5_fa_nf_pairs`; dossier sec. 2.5) | Cys 33 + sugar 100 mM, 0.5 M phosphate pH 5, 145 C, 100 mL in a 200 mL autoclave: NF / MFT 2754 (ribose), 7.1 (glucose), 23.9 (rhamnose); furfural / FFT 5.58 / 2.46 / 7.75. NF and furfural method not stated (MFT, FFT by SIDA) | whether glucose MFT passes through the intermediate that ribose MFT does (D1) | DEV; ribose NF / furfural are fit rows; the glucose cells run under R0 only and are spent (ruling 2) |
| Hofmann 1998 Table 9 (`hofmann1998_t9_m2p_glucose_cys`) | 1-mercapto-2-propanone from glucose + Cys 11.4 / 59.6 / 26.5 ug at pH 3 / 5 / 7 | a C3 marker in the same pots, maximum at pH 5 | DEV; **seen** (`thiol_submodel.yml`) |
| de Roos 2005 Table II (`deroos2005_table2_cysamadori_glucose_aa`) | glucose + Cys + Gly or Pro, 90 C: D-fructose-L-cysteine forms readily; Fru-Cys / Fru-Gly 0.9 - 1.1 (15 % water), 4 - 5 (79 % water) | the Cys-hexose Amadori compound exists in situ; amine competition for glucose | DEV (Givaudan); relative RI areas |
| de Roos 2005 Table I (F23#11; `deroos2005_table1_thiol_yields_cys_carbonyl`) | 1 : 1 Cys + sugar, lactate pH 5, 37 % water, boiling 1.5 h: MFT ribose 43.5, xylose 36, glucose < 0.05 mg/mol Cys (the paper's reporting limit) | at ~100 C the hexose route is shut although its Amadori compound forms (Table II) | DEV |
| **Du 2024b** (CAAS, Wang group; `du2024b_extraction.md`; new here) | Cys + glucose 200 / 200 mM, 0.2 M phosphate pH 6.5, 145 C (no-aldehyde control). The Cys-Amadori compound peaks at ~124 mM (~62 % of the cysteine) at 10 min and falls to ~5.5 g/L by 30 min. Cys + Cys-Amadori exceeds the charge by 40 % (a mass-balance defect). No MFT or NF in the control | the Cys-hexose Amadori pool forms within minutes at Hofmann's temperature, and decays at >= 0.06 - 0.16 /min | DEV; a screen signature (S5, sec. 5.5) |
| **Li 2011** (`li2011_extraction.md`; held, cited here for the first time) | glucose + Cys ~0.2 M each, pH 7 phosphate, 100 / 110 / 120 / 130 C: free-cysteine first-order loss 0.0011 / 0.0032 / 0.0081 / 0.0203 /min, Ea 121 kJ/mol | the net formation of the Cys-hexose pool. [dry] At 145 C its Arrhenius line gives 0.074 /min; Du 2024b's peak implies >= 0.097 /min of free-Cys loss, x1.3 apart, inside Du's mass-balance defect | DEV; **seen** (B59 hold-out diagnosis); a band source |
| Zhang 2026 (`zhang2026_fig3_gsh_glucose_nonvolatiles`, `zhang2026_tables6_s8_volatiles`) | GSH + glucose 60 / 60 mM, 0.2 M phosphate pH 6.5, 130 - 150 C: Cys-, GSH-, Gly-Amadori by LC-MS/MS (figure-only); MFT only at 90 / 120 min (0.33 - 2.32 ng/g, RF 1) | the Cys-hexose Amadori pool in a peptide pot; MFT's late onset | DEV (BTBU); **seen** |
| Wang 2020 (BTBU; `wang2020_extraction.md`; held, cited here for the first time) | GSH + glucose 60 / 60 mM, 0.2 M phosphate pH 6.5, 140 C, 0.16 - 4 h: the GSH-Amadori compound and Cys-Gly-sugar thiazolidines over time; no MFT | the peptide-hexose Amadori pool (part A) | DEV |
| Hofmann 1999 (F06#4, F07#20) | xylose osones peak at 10 min, glucose at 40 min | the hexose supply's time scale | DEV |
| Brands 2001 (F14#0 - #1, F02#13, F05#0) | 3 % caseinate + glucose, 120 C: net glucose loss identical with and without casein; formic x1.4, acetic x2.6 with casein | protein-bound lysine changes products, not sugar loss | DEV |
| Xia 2022 (F02#4 - #5) | Glc + Gly / diglycine / triglycine, separate pots: TriG-ARP forms 6.9x (100 C), 5.0x (130 C) faster | peptide amine as its own class | DEV (Ho network) |
| **Cannon 2023** (Rutgers, Ho; `cannon2023_extraction.md`; new here) | glucose + Cys.HCl 0.1 M each, water, label pH 7.4 (final 4.4), 1 h: MFT 0.43, FFT 0.33, MP 3.97 mg/mol Cys at 160 C; all three nd at 130 C (no LOD); octanol IS on FID, no recovery; **no NF**. Arg arms | part A (Arg vs Cys at fixed glucose); the hexose route shut at 130 C and open at 160 C in one study | DEV |
| Du 2025 (BTBU) | GSH + glucose 200 / 200 mM, 0.2 M phosphate pH 6.5, 145 C: MFT 135 ug/L at 120 min (5.9e-6 mol/mol GSH); CAMOLA on five sulfur volatiles, **not on MFT** | MFT level in a peptide pot (part A, weak); its lipid arms are M5's | DEV |

### 2.5 Mass balance

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Brands 2001 (F05#1, F02#12) | fructose + casein: lysine lost without Heyns or furosine | the protein amine's fate is not the free-amine pattern | DEV |
| Cao 2017 (F04#17, F04#18) | the Cys-bearing pool conserved over 5 h; glycine at 16 - 24x Cys drains TTCA to ~0.03 M | amine competition in a cysteine pot | DEV (BTBU); seen (M1c) |
| Zhang 2025b (F04#19) | lysine at 16x Cys moves TTCA and Cys-Amadori | the same with a free basic amine | DEV; seen (M1c) |
| Kanzler 2026 Fig. 1F | amino acids left after 120 min with MGO: Ala 100, Lys 98, Glu 90, Phe 63, Arg 47.5, Gln 15 % | which amine classes are consumed by a dicarbonyl and which only catalyse its loss | DEV |
| Yu 2017 (Jiangnan; `yu2017_extraction.md`; new here) | glucose + Lys / Arg / Gly / Pro 40 mM, 120 C, 0.5 h: glucose -66 to -85 %, amino acid -3 to -30 % | amine catalysis of sugar loss without matching consumption (part A) | DEV |

### 2.6 One factor varied inside one study

| source (record) | factor -> result | bears on | split |
|---|---|---|---|
| Hofmann 1998 Table 1 (F23#21, F25#10; `hofmann1998_t1_cys_carbohydrates`) | ribose / xylose / fructose / glucose / rhamnose + Cys, pH 5: MFT 19.8 / 14.3 / 2.5 / 1.9 / 0.8 ug; FFT 12.1 / 9.6 / 3.2 / 2.8 / 0.8 | ribose / glucose MFT 10.4, FFT 4.3 (D4) | DEV; **seen** |
| Hofmann 1998 Table 2 (F23#18, F25#8; dossier sec. 2.2) | 0.5 M phosphate pH 3 / 5 / 7: glucose MFT 0.3 / 1.9 / 0.4, FFT 0.7 / 2.8 / 0.6; ribose MFT 55.3 / 19.8 / 2.5, FFT 22.9 / 12.1 / 1.2 | the hexose route's pH shape against the pentose one (D2) | DEV; glucose pH 3 / 7 are Maillard-path hold-out bundles, **run under R0 only; spent** (ruling 2), never fitted (ruling 3) |
| **Ames 2001** (Reading; `ames2001_extraction.md`; held, cited here for the first time) | Cys + glucose 0.044 mol/kg in wheat starch, 18 % moisture, twin-screw extrusion, label pH 5.5 / 6.5 / 7.5 x die 120 / 150 / 180 C. MFT (ng/10 g): 180 C nd / 7 / 31; 150 C 6 / 3 / 37; 120 C nd / 6 / 15 at pH 5.5 / 6.5 / 7.5: **rising with pH**. Single runs; measured product temperature departs from the die setting by up to +15 % | **against** the fall R1 / R2 / R5 predict on the 5 -> 7 limb, if it transfers to water. It probably does not (low water, starch, die-temperature distortion; P(transfer) ~0.2) | DEV (v2 E8 lab: Reading); not simulable (no extrusion history or water basis in the network) |
| Meynier 1995 Table 3 (F25#9, F13#5; CF-001 conditions) | Cys 50 + ribose 35 mM, 0.2 M pyrophosphate, held (cooled drift <= 0.2), 140 C, 1 h, pH 4.5 -> 6.5: FFT 8.5 -> 1.4 (x6.1), furfural 3.9 -> tr (> 39x), MFT 15.2 -> tr | E8: FFT falls far less than furfural (D6) | DEV(E8) |
| Hofmann 1998 Table 2 (F25#8) | ribose FFT x19 from pH 3 to 7, highest at pH 3 | E8 | DEV(E8) |
| Cerny 2007 (F25#7, F13#7) | FFT x2.3 down from pH 4 to 6; zero at pH 7 | E8, third lab | v1 DEV; v2 excluded (calibration input) |
| **Yeo & Shibamoto 1991** (UC Davis; `yeo1991_extraction.md`; new here) | Cys + glucose 1.0 M each, water, label pH 2 / 5 / 7 / 9, microwave 15 min (temperature not measured): FFT area ratio 0.18 (pH 5), 0.22 (pH 7), nd at 2 and 9; **no MFT, no NF** | part E only (FFT vs pH, hexose); not D2 | DEV |
| Schirle-Keller 1992 (Minnesota) | glucose 0.5 M + Cys 0.1 M, unbuffered, 80 - 150 C: net furfural rises ~30x from 100 to 120 C; printed apparent Ea 147 kJ/mol (refit 156); no pH, no sulfur volatiles | `r_glc_fur`'s temperature law (part E; deferred) | DEV |
| Mottram 2002 (F23#20, F14#8) | buffered vs unbuffered: MFT x2.2 - 2.6, FFT unchanged; ribose-only furfural unchanged | the branches respond differently (shared with M5) | DEV |
| Blank 1997 (F13#11, F14#10) | HDMF from xylose needs glycine or alanine (< 0.01 vs 2.6 ug/mmol) | amine catalysis at 90 C | DEV; calibration input |
| **Hayashi & Namiki 1986** (`hayashi1986_extraction.md`; new here) | glucose + beta-Ala 2 M each, water, initial pH 3.51 / 6.42 / 9.29, ~100 C: C3 fragment ~0 / 0.013 / 0.117 mM peak; C3 initial rate 0.41 dex/pH (6.4 -> 9.3), acid limb >= 0.28 dex/pH (censored); C2 0.23 dex/pH | R3's fragment-supply pH law (a band, sec. 5.2) | DEV |
| **Yu 2017** (as sec. 2.5) | 0.2 M citrate-phosphate pH 5.0 / 6.0 / 7.0 / 7.4 / 8.0, 120 C, 0.5 h: MGO rises 0.50 - 0.76 dex/pH (5 -> 7), 0.18 - 0.36 (7 -> 8); at fixed pH Lys gives 2 - 7x Gly's MGO, Pro 0.2 - 0.8x | R3's C3 supply vs pH on D2's limb (a band); within-class amine spread as large as D7's between-class contrasts | DEV |
| **Luo 2022** (Ho network; `luo2022_extraction.md`; new here) | glucose + glycine + diglycine (or triglycine), 20 mM each **in one pot**, initial pH 3 - 9, 80 - 160 C, 1 - 2 h: diGly-ARP / Gly-ARP and triGly-ARP / Gly-ARP 0.43 - 2.33 at 80 - 100 C, pH 3 - 7 (geometric centre ~1); 0.08 - 0.17 at pH 9; 0.02 - 1.15 at 130 - 160 C | **part A against Xia 2022's peptide advantage** in competition at fixed sugar: equal per alpha-amine (A-R3) at the low-T, acid-neutral corner | DEV (Ho network, as Xia) |
| Tang 2024 (F05#16 - #17) | soy vs pea isolate, wet 85 C: lysine -22 vs -8 %; the order flips dry | protein identity as a class constant | DEV (Wageningen) |
| Zhou 2023c, Feng 2022 (F27#24 - #26) | GSH a weaker MFT source than Cys; the ribose-GSH ARP makes MFT itself | peptide vs Cys class | DEV; Feng a calibration input |
| Zhang & Ho 1991b (Rutgers) | glucose + Cys vs glucose + GSH, 50 / 50 mM, water, label pH 7.5, 180 C, 1 h: total volatiles 969 vs 404 mg/mol; **no MFT or FFT detected in either pot** (no LOD) | part A (Cys vs peptide at fixed sugar); a non-detect of glucose MFT at pH 7.5, 180 C | DEV |
| Tai & Ho 1997, 1998 (Rutgers; held, cited here for the first time) | 1997: Cys + glucose equimolar, water, pH 6, 160 C, 2 h: sulfur volatiles, **no MFT reported** (no LOD). 1998: GSH + glucose, initial pH 3 / 6 / 8, 160 C, 2 h: FFT 9 - 15 mg/mol, nearly flat in pH; **no MFT** | part A and part E; glucose MFT non-detects (no LOD) | DEV |
| Martins 2010 (Unilever) | 32-run D-optimal design, 100 C, 0.2 M pyrophosphate pH 6: xylose, glucose, 5-ketogluconate, Cys, Gly, Glu; MFT, FFT, NF responses as regression signs only | part A competition (Gly vs Cys at fixed sugars), weakly; no per-run data | DEV |
| Bolton 1993 (CF-090) | + 6.3 mM xylose leaves MFT unchanged in the glucose matrix | thiamine dominance in hexose foods | DEV |
| Lee 2022 (F06#16), Kocadagli 2016 (F06#15, F07#18) | amino acid raises osones 4 - 10x; NaCl moves the split | amine catalysis | **v1 TEST**: listed, not used |
| Xu 2007 (`heldout_xu2007_tab1.yml`) | pH 5.6 / 7.0 / 8.0 control arms | E8 | keyed held-out (not in a frozen split): item 7 only |

### 2.7 Fitted time courses and fits

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Zhang 2026 Table 2 (`zhang2026_table2_rate_constants`) | 19 k at 130 / 140 / 150 C; Cys + glucose -> Cys-Amadori Ea ~101, Cys-Amadori -> alpha-dicarbonyls 85.7 kJ/mol (both already B59b's transferred barriers); units not printed, no CI | barriers for part G | DEV (BTBU) |
| Martins 2005, 2005b (F06#12 - #13, F07#15 - #16) | k4 +0.35, k7 +0.60 dec / pH; Ea 97.1 / 107.3 | the hexose split | DEV; calibration input |
| Brands 2002c, Nguyen 2016, Wolf 1977a (F05#3 - #8; CF-072, CF-096, CF-100) | glycation 0.39 dex / pH; Ea 97 - 138 kJ/mol; lysine return disputed | the protein-amine class constants | DEV |
| engine records | B8, B9, B47, B51.1, B20, B21 (`k_ama_g`, `k_g_go`), B1 (`k_mgo_mel`), EX-012, EX-013, red-team probe | sec. 1 table | engine records |

### 2.8 Left out, with reasons

- The thiamine route (F27, E18), except Bolton's bound: a separate MFT source with its own case.
- Acrylamide amine competition (Claeys 2005, Rydberg 2003; F29#3, #21; De Vleeschouwer 2009): E9, the acrylamide lane.
- Pyrazine amine identity (Leahy Lys >> Asn, F19#10; Van Lancker dipeptides): M6.
- Strecker aldehydes from glucose pots: M3. Lipid with glucose (Du 2025's fatty-acid arms, Zhao 2019b Table 1, Yang
  2015, Negroni 2000): M5.
- Goncuoglu 2026 PPI extrusion (F29#24, v2 TEST E9) and Hamzalioglu 2026 (F05#10 - #11, v1 TEST): protein amine in
  acrylamide / milk chemistry; held out. R5 uses B21's shipped `k_ama_g` (fitted on Hamzalioglu's rates) at its shipped
  value only; no Hamzalioglu value is read here.
- Hofmann 1997 (glucose / rhamnose + Cys, the Table 1 pots): FD factors only, no amounts. Hofmann 1996b (Weurman):
  declared duplicate of Hofmann 1998's tables (`pdf_dossier_aliases.yml`).
- Schieberle 2005 Table 3 (fed C3 + C3 -> HDMF): two of three cells repeat Hofmann 1998 cells exactly; the primary is
  Frauendorfer & Schieberle (not on disk) and HDMF is not an M4 species. Hofmann 1998 Table 6 is used instead.
- Madruga 1995 (beef psoas, pH 4 - 5.6, MFT rises in acid): a meat matrix whose sugar is mostly IMP-derived ribose, not
  a hexose route. W5 (recipe) evidence, not M4's.
- Liu 2021 (FFT from glucose + Cys with ribose); US 9,005,689 (glucose + Glu + Cys, FFT only, no method, a patent):
  part E, deferred.
- El-Ghorab 2010: no volatile analysis is described; its MFT statement is a qualitative discussion sentence. Nothing to
  list.
- Fanun 2001: not M4 (a Vauthey sibling, passed to H0's question); it prints no amount.
- Haleva-Toledo 1999 (70 C, furaneol + Cys, no MFT), Zyzak 1995 and Huyghues-Despointes 1996 (physiological or
  pyrolysis conditions, no thiol at process temperature): outside the criterion.

**Counts by kind (sources).** Labelling 11 (one gap); fed intermediates 13; trapping 6 (two held out); measured
intermediates 13; mass balance 5; one factor within a study 21 rows (three held out); fitted 4 plus the engine records.

### 2.9 The citation sweeps (sec. 7 item 2)

**The criterion each sweep was keyed to.** The question asked of each work: does it hold
- (a) NF and MFT measured in one hexose + cysteine (or H2S) pot (D1);
- (b) hexose + Cys MFT at >= 2 pH in one study (D2);
- (c) a positional glucose label on MFT, or a fed hexose intermediate (1-DG, the Cys-hexose Amadori compound) giving
  MFT;
- (d) an amine-class contrast at fixed sugar (D7);
- (e) a constant for a band this case needs (the Cys-hexose Amadori steps, 1-DG's exits, MGO's loss law)?

**Added 2026-11-06 with R5:**
- (f) an O2 (air vs inert) contrast on hexose + Cys MFT, or on hexose C5 dicarbonyls at >= 100 C;
- (g) a fed glucosone giving NF, MFT or a C5 dicarbonyl;
- (h) a positional glucose label on NF or 1-DPO at >= 100 C.

Screening used titles and abstracts only (blinding holds).

**Sweep 1 (2026-11-05), criteria (a) - (e).** Run with `scripts/tools/openalex.py sweep` (from `origin/openalex-tool`,
PR #70) through the wrapper, ~1 request per 1.1 s, anonymous.
- Seeds, forward and backward: Hofmann 1998 (10.1021/jf9705983), Du 2023 (10.1016/j.foodres.2023.113337), Zhao 2019b
  (10.1016/j.foodchem.2018.07.106), de Roos 2005 (10.1021/bk-2005-0905.ch009), Zhang 2026
  (10.1016/j.foodchem.2026.148681), Tressl 1994 (10.1021/bk-1994-0564.ch018), Voigt 2009 (10.1021/jf900459x).
- Seeds, forward only: Meynier 1995 (10.1016/0308-8146(95)93282-v), Xia 2022 (10.1021/acs.jafc.2c06639), Kanzler 2026
  (10.1016/j.foodchem.2026.151175), Davidek 2006a (10.1021/jf060667q), Deng 2021 (10.1021/acs.jafc.1c02023).
- Result: 634 distinct works, 123 already held, **511 unseen**. Decisions:

| decision | n |
|---|---|
| food or matrix product study, no defined hexose + amine pot | 194 |
| non-thiol chemistry outside the criterion (colour, antioxidant, acrylamide, HMF, pyrazines, dicarbonyl surveys, ARP synthesis) | 72 |
| read individually: no hexose + Cys / amine pot bearing on (a) - (e) | 74 |
| lipid interaction (M5) | 48 |
| biomedical or other field | 43 |
| review, book or chapter without new pot data | 40 |
| pentose-only system (M1 / M2) | 22 |
| computed (DFT, ReaxFF): excluded by the no-computed-barrier rule | 3 |
| held by file name, not DOI (Hofmann 1996b, a declared duplicate) | 1 |
| not M4 (Fanun 2001), passed to H0's question | 1 |
| candidates | 13 |

Notable exclusions:
- Yaylayan & Keyhani 2000, 2001: positional [13C]glucose with Ala / Met, but pyrolysis GC-MS, no thiol at process
  temperature.
- Goto 1961 / Hayami 1961: 14C acetol, alkaline sugar chemistry.
- Shibamoto & Russell 1976, 1977: glucose-H2S-NH3, qualitative volatiles.
- Lee 2010: GSH + glucose CAMOLA on two thiophene aldehydes, not MFT; FD factors.
- Song 2019: [13C6]glucose + Cys + 2-nonenal, MFT not named; M5.
- Yang 2015: Cys + glucose with fat; M5.
- Cerny 2008, Mottram 1994: reviews of work already held.

**Sweep 1's candidates, read 2026-11-05 under roles fixed before reading (PR #76; all first reader).**

| candidate | role fixed before reading | outcome |
|---|---|---|
| Yeo & Shibamoto 1991 (10.1021/jf00002a029) | decisive for D2 (second lab) only if MFT at >= 3 pH | **NOT MET**: no MFT; FFT at 4 pH (part E evidence, sec. 2.6) |
| El-Ghorab 2010 (no DOI) | as Yeo | **NOT MET**: no volatile analysis (sec. 2.8) |
| Cannon & Ho 2023 (10.31665/jfb.2023.18338) | decisive for D1 (second lab) only if NF and MFT in one pot with recovery | **NOT MET**: no NF, no recovery; part A evidence (sec. 2.4) |
| Du 2024b (10.1016/j.lwt.2024.115741) | band source for the Cys-hexose Amadori steps, by amendment; never decisive | met with a defect (mass balance +40 %); **used as screen signature S5 instead of a band** (sec. 5.5), so the bands stay disjoint from the screen |
| Schieberle, Fischer & Hofmann 2003 (Weurman 2002, no DOI) | evidence for D3 only | **not obtained**: the file read was Schieberle 2005 (10.1196/annals.1333.029), a review (sec. 2.1). What the case loses: sec. 4.5 |
| Tressl & Rewicki 1999 (10.1007/978-1-4615-4693-1_26) | decisive for the which-carbon gap only if MFT isotopologues from a positional label | **NOT MET**: no MFT datum; Fig. 14 conflicts with Tressl 1994 (sec. 2.1) |
| Glomb, Gobert & Voigt 2010 (10.1021/bk-2010-1042.ch004) | band source for 1-DG's exits and `r_odg_dpo`'s ceiling; never decisive | met qualitatively (no floor); **names rival R5** (sec. 3) |
| Hayashi & Namiki 1986 (10.1080/00021369.1986.10867692) | band source for R3's fragment supply pH law | met; band (sec. 5.2) |
| Pfeifer, Haase & Kroh 2013 (10.1021/jf302959k) | band source for MGO self-loss, second lab to be confirmed | one temperature only; **same laboratory lineage**, no second lab |
| Kanzler, Wustrack & Rohn 2021 (10.1021/acs.jafc.1c04594) | band source (MGO + Lys); same lab | met at 100 C only; **k2 x3.2 Kanzler 2026's**: the MGO-alone Ea basis falls to 47 - 66 kJ/mol |
| Yu, Xu & Yu 2017 (10.1002/jsfa.8158) | part A evidence; band for MGO supply vs pH | met in part; band (sec. 5.2) |
| Luo, Xu, Li & Ho 2022 (10.31665/jfb.2022.18332) | part A evidence (D7) | met; **points the other way from Xia 2022** (sec. 2.6) |
| Luo et al. 2023 (10.1021/acs.jafc.3c06250) | band source for peptide-ARP decay vs T | met; apparent Ea 34 - 56, below the band; Cys-return low end 0.54 (sec. 5.2) |

No decisive role was met. D1 and D2 keep one laboratory (DFA), and the which-carbon gap stays open.

**Sweep 2 (2026-11-06), keyed to criteria (f) - (h).** This sweep was run because R5 changed the decisive criterion's
wording after sweep 1.
- Seeds, forward and backward: Glomb 2010 (10.1021/bk-2010-1042.ch004), Gobert 2009 (10.1021/jf9019085), Tressl &
  Rewicki 1999 (10.1007/978-1-4615-4693-1_26).
- Seeds, forward only: Cannon 2023 (10.31665/jfb.2023.18338), Du 2024b (10.1016/j.lwt.2024.115741), Smuda 2010
  (10.1021/jf100334r), Voigt 2010 (10.1021/jf100140h).
- Result: 375 records, 353 distinct works, 40 already held by DOI, **313 unseen**. Every unseen work was screened on
  title and abstract: the 83 with three or more criterion terms individually, the rest by title. Overlap with sweep 1's
  unseen list is not recorded; its table was not committed.

| decision | n |
|---|---|
| other Maillard chemistry outside the criterion (pyrazines, pyrrolizines, AGEs, ARP synthesis, colour, antioxidants) | 107 |
| food or matrix product study, no defined hexose + amine pot | 78 |
| biomedical or other field (in-vivo dicarbonyls, dialysis fluids, lens proteins) | 73 |
| dicarbonyl or fragment chemistry with no thiol, no O2 contrast at >= 100 C and no fed C6 dicarbonyl | 34 |
| review, book or chapter without new pot data | 13 |
| sulfur chemistry without a hexose + Cys MFT / NF pot (pentose ARPs, Met, diacetyl, thiamine) | 3 |
| **candidates** (below) | 5 |

Notable exclusions:
- Yaylayan & Keyhani 2000 (10.1021/jf000004n): pyrolysis, as in sweep 1.
- Werkhoff 1990 (10.1021/jf00093a041): the thiamine route.
- Hayashi 1986 companion (10.1080/00021369.1986.10867691): C3 fragments, for R3, which is a declared-failing control.
- Glomb's milieu-isomer NMR papers (10.1016/j.carres.2012.10.010, 10.1021/jf305302e): structure, not route.
- Hollnagel & Kroh 1998 (10.1007/s002170050294): fragments from sugars.
- The Glomb-group fragmentation review (10.1021/bk-2016-1237.ch009): a review.

**Sweep 2's candidates.** Each is unseen. Its role is fixed now, before reading. P = the chance it holds the named
measurement.

| candidate (from the OpenAlex record) | DOI | role, fixed now | the measurement that would count | P |
|---|---|---|---|---|
| "Analysis of volatiles from the thermal decomposition of Amadori rearrangement products in the cysteine-glucose Maillard reaction and density functional theory ..." (Food Res. Int. 2024; title only; cites Du 2024b) | 10.1016/j.foodres.2024.114454 | **decisive for D1 in its fed form only if** NF and MFT are both quantified from the fed Cys-glucose ARP with stated recovery, and the laboratory is not one whose pots are M2's TEST (BTBU, CAAS); otherwise evidence for D5 and arm 3. Any computed (DFT) part is excluded | NF and MFT from a fed Cys-glucose ARP, charges and T stated | MFT reported ~0.4; decisive ~0.1 |
| "Browning Potential of C6-alpha-Dicarbonyl Compounds under Maillard Conditions" (JAFC 2017) | 10.1021/acs.jafc.6b04512 | band source for R5's glucosone -> C5 step and for `r_odg_dpo`'s ceiling at process temperature (fed glucosone, 1-DG, 3-DG + Ala, 130 C, pH 5 / 8); never decisive (no thiol) | C5 dicarbonyls (1-DPO, 3-DPO) per fed C6 dicarbonyl | ~0.35 |
| "Mechanism of pyrazines and thioethers formation promoted by high oxygen concentration in the methionine-glucose Maillard reaction system" (JSFA 2024) | 10.1002/jsfa.14094 | band source for the O2 law of hexose dicarbonyl supply (R5) at process temperature; never decisive (Met, no MFT) | a dicarbonyl or osone level at >= 2 O2 levels | ~0.3 |
| "Formation of Reactive Intermediates, Color, and Antioxidant Activity in the Maillard Reaction of Maltose in Comparison to d-Glucose" (JAFC 2017) | 10.1021/acs.jafc.7b04105 | band source: glucose + Ala C5 dicarbonyl and glucosone standing levels at 130 C, pH 5; never decisive | 1-DPO / 3-DPO / glucosone vs time from glucose | ~0.3 |
| "Dicarbonyl-mediated co-formation of aroma compounds and AGEs across heating times" (Food Chem. 2025; title only; cites Du 2024b) | 10.1016/j.foodchem.2025.146741 | evidence only, if it is a hexose + Cys pot with dicarbonyls and MFT over time; never decisive (likely CAAS) | MFT and C5 dicarbonyls over time in one hexose + Cys pot | ~0.1 |

The chance that any sweep 2 candidate supplies a decisive unit is ~0.1. None sets a band in this case until it is read
and entered by an amendment registered before any M4 engine result exists.

**Roles confirmed before reading (2026-11-06; committed before either paper or its SI was opened).** The owner has
obtained two candidates and their supplementary material. No M4 engine result exists. Their roles stand as already
declared and are not changed after reading:
- **10.1016/j.foodres.2024.114454** (the fed Cys-glucose ARP volatiles paper; sweep 2's first candidate):
  - **decisive for D1 in its fed form only if** NF and MFT are both quantified from the fed Cys-glucose ARP with
    stated recovery, and the laboratory is not one whose pots are M2's TEST (BTBU, CAAS);
  - **otherwise evidence** for D5 and the arm-3-type question (MFT per mol of fed Amadori; R4' vs the claim);
  - any computed (DFT) part is excluded;
  - P(MFT reported) ~0.4, P(decisive) ~0.1.

  Under ruling 2 (falsify, never credit) even a met decisive role could only falsify a constant-free prediction. It
  could not credit a rival on the basis of a paper read after the case named it.
- **Wang Y. et al. 2019, 10.1016/j.foodres.2019.05.001** (M0a's candidate, carried into part A; sec. 1):
  - **decisive for part A (HMF) only if** an Asn + hexose pot carries HMF at >= 2 pH with the lipid held fixed;
  - otherwise evidence;
  - part A is structure only (FM7), so a met role feeds a later part-A case, not this one;
  - P ~0.5.

The SI of each is read under the same role. A dossier for each records whether the role is MET.

**Outcomes (read 2026-11-06 after the confirmation above; first reader).**

| paper | dossier | role | outcome |
|---|---|---|---|
| Lei et al. 2024, Food Res. Int. 188:114454 (Dalian Polytechnic University) | `lei2024_extraction.md` | decisive for D1 (fed form) only if NF and MFT from the fed Cys-glucose ARP with recovery, not a BTBU / CAAS lab; else evidence for D5 / arm 3 | **NOT MET.** Purified "ARP" (identity not secured: NMR anomalies), 0.01 mmol, **dry**, air, micro-chamber at 120 - 180 C for 120 min. Semi-quantified against one internal standard (RF 1), no recovery. MFT only at 150 / 160 C (4.00 ± 1.11 / 2.04 ± 0.03 ug/g), identified by RI / MS only, with RI 979 vs the cited LRI 869. **No norfuraneol.** No consumed-ARP denominator, no matched glucose + Cys pot. The SI is DFT rate constants only (excluded). Fallback: weak evidence for D5 (a fed pool gives sulfur volatiles and a library-matched MFT). It does not separate R4' from R1 / R5 |
| Wang Y. et al. 2019, Food Res. Int. 123:403 (Nanchang University) | `wang2019_extraction.md` | part A: decisive for HMF only if an Asn + hexose pot carries HMF at >= 2 pH with lipid fixed | **MET IN FORM, CONFOUNDED.** The lipid-free Asn + glucose control at 180 C, 30 min gives HMF ~12 (water, initial pH 4.6) vs ~9 nmol/mL (0.4 M acetate, pH 6.0) (Fig. 2C, read by eye). But buffer and pH change together, and the cooled pH converges (5.7 vs 5.9; SI Table S7). Same pots: 3-DG 65.6 vs 24.6, MGO 65.7 vs 322 nmol/mL (Table 2). It goes to the later part-A case with a declared acetate term and an explicit pH scale; part A is structure only here |

**Does either change the decision? No.**
- Neither supplies a credited unit for part G: Lei misses NF and recovery, and under ruling 2 nothing read after the
  freeze-ready text could credit a rival anyway.
- Neither touches R1 vs R5 (no label, no O2 contrast), R4' (no arm-3 denominator) or R0's falsified route.
- Wang 2019 is part-A evidence, and part A is an engine-structure requirement only (FM7).
- So the status stays "READY TO FREEZE as development grade, decider W6, pending the reduced-model power check". The
  plan stays: build only the W3 structure (part G's species and the part-A classes, at parity, off by default) and the
  reduced forced-input mini-model.
- Neither paper becomes a band source; neither enters a ledger or a fit table.

**Held dossiers (2026-11-06 screen).** Open question 7 asked for a screen of every held dossier that names glucose
with a thiol or a deoxyosone. That screen found 199 such files; 135 were not cited in the 2026-11-05 text. Each was
screened on its header and on keyword counts by one reader (not two blinded readers, as M1c did). Those that bear on
the criterion now sit in sec. 2 or 2.8 with their reason: Ames 2001, Li 2011, Tai 1997 / 1998, Wang 2020, US 9,005,689,
Haleva-Toledo 1999, Zyzak 1995 and Huyghues-Despointes 1996. The rest are pentose-only pots, thiamine pots, meat and
other food matrices, patents on heme systems, trunk-only dicarbonyl kinetics (Goncuoglu, Gokmen, Yu 2018 / 2020,
Zhang 2020), downstream fed furanones (Zheng 1997, US 3,904,655, Munch 1997 / 1998), and reviews. Ames 2001 is the one
the case had missed that speaks against a rival (sec. 2.6).

### 2.10 Spent data (sec. 7 item 2)

**The rule (FM2, adopted 2026-11-06).** Spent or seen data may **falsify** a rival whose prediction there is
constant-free (structural: a label share a route must leave, a species a structure lacks). They never **credit** a
rival. The rule applies to R0's route exactly as to the new rivals.

**Spent: these data nominated a rival, and are development evidence for it.**

| data | spent by | spent for | consequence |
|---|---|---|---|
| Hofmann 1998 Table 1 (all arms) and the T1 design | atlas E1; probes C and E; the red-team probe; M0 CK-41 (MF-12's size); B8 - B51.1 level rows | R1, R4 (the claim and the band-now fix were nominated on it) | D4 is side evidence only; S3 is a screen signature |
| **Du 2023 / Zhao 2019b label shares (D3; F23#7 / #8)** | **the data that nominated E1** (ENGINE_VS_LITERATURE l. 259 / 298; corrected by FM2) | the claim (the hexose route is an intact C5) | **falsify-only:** D3 falsifies C2 + C3 as the hexose MFT route (R0's route and R3; constant-free), credits nothing. Label arithmetic only; no engine run on these pots (M2's TEST) |
| **Hofmann 1998 Table 5 glucose NF / MFT (D1)** | seen through R0's no-refit probe | -- | **falsify-only** development unit: reported, never credited |
| **Hofmann 1998 Table 2 glucose MFT at pH 3 / 7 (D2)** | MF-12 observations (M0); R0's probe | -- | **falsify-only** development unit; never fitted (ruling 3) |
| Bolton 1993 / 1994 | core panel; M0 MF-15 | -- (no M4 rival was nominated on it) | S4 only |
| Glomb 2010 Table I and its label statements | the 2026-11-05 dossier, which nominated R5 | R5 | never decisive (it is not a thiol pot anyway) |
| Steinhagen 2021 (glucose alone, 121 C) | M0 MF-23 (read for sec. 2.11) | R5's oxidative piece | sign evidence only |
| hofmann2000 Fig. 2, Tables 3 - 4 (glucose + Phe PA) | M0 CK-45 (planned), MF-30 | R2's 1-DG + amino-acid Strecker piece | sign evidence only (sec. 2.11) |
| Kanzler 2021 / 2026, Pfeifer 2013, Bork 2022 | CF-122, M0 MF-24, CK-46 / CK-47 (planned) | the MGO loss law (a band, not a rival) | band source only, never scored |
| Hofmann 1998 Tables 6, 8 | the 2026-11-05 R3 correction | R3 | fixed measured constants for R3, never scored |
| **Yu 2017** (glucose + Lys / Arg / Gly / Pro, MGO vs pH, 120 C) | this case (R3's fragment pH band; part A's within-class spread; sec. 2.11's MF-34 reading) | R3, part A | **SEEN / spent: never credited as decisive anywhere** (here, in M6 or in M0a). Band and evidence only |
| Schibilsky 2019 (glucose + Ala DMHF, two pH arms) and M0a's other spent data | M0a (CK-06, sec. 2.8), now folded into part A | part A's pH module | spent for part A |

**Seen but not spent.** These have engine numbers on record but nominated no rival here.
- Li 2011, Deng 2021, Martins, Davidek 2006a, Zhang 2026: seen in other records; band sources here, never scored.

**Spent by the screen once it runs:** de Roos 2005 Tables I and II, Hofmann 1998 Table 1, Bolton, Du 2024b.

### 2.11 M0's groups: what an M4 rival predicts for their members on seen data

The protocol's best evidence for a mechanism is a shared piece that predicts a misfit beyond the one it came from.
These are the members of the two groups M0 assigns to M4, checked against each rival's new pieces. Directions are
structural or [dry]. No engine was run.

**G4 (the supply reservoir: MF-01, -02, -03, -04, -05, -15).**
- MF-01 to MF-05 are pentose + Cys pots (Schieberle 2000, Cao 2017, Zhang 2025b, Zhai 2019 / 2021, Zhao 2017b). Part G's
  new steps need a hexose, so **every M4 rival is inert there and predicts nothing**.
- The shared piece is structural only. The Cys-hexose Amadori pool (R1, R2, R5) is the hexose member of the class
  whose pentose member is M1c's reservoir (`r_cama_form`, `CYS_AMADORI_ROUTE`).
  - Rule declared here: **M4's bands never come from M1c's fitted reservoir constants.** Those are fitted on G4's
    pots and would import their compensations (CF-107).
  - If M1c and M4 both ship a reservoir, the joint refit judges the class constant once.
- MF-15 (thiamine pots), measured basis: Bolton 1994's glucose bouillon is the one member with a hexose. There
  R1 / R2 / R5 sequester Cys in the hexose Amadori pool and add a small Cys-derived MFT route. The engine's thiamine MFT
  has no cysteine term (MF-15 candidate (a)), so the +1.30 dex miss barely moves: no discriminating prediction. Bolton
  is S4 (a screen) in any case.
- **Verdict: G4 gives M4 no protocol evidence.** It gives one dependency (the class) and one rule (no bands from
  M1c's fits).

**G8 (hexose dicarbonyl / osone supply: MF-12, -23, -30, -34).**

| member | rival with a piece that reaches it | prediction on the seen data | reading |
|---|---|---|---|
| MF-12 (hexose MFT ~0) | all part-G rivals | M4's own misfit (CK-41 done; spent) | not evidence beyond itself |
| MF-23 (glucose alone, 121 C: glucosone -1.80, glyoxal -1.54 dex) | R5 | R5 as specified acts on the Amadori compound, so in an amine-free pot it **predicts nothing**. If its beta-cleavage step sits on the one glucosone pool (as it must on the network), it diverts glucosone from `r_g_go` -> glyoxal and makes the glyoxal miss **worse**; the size is set by the step's band (0 to a ceiling) | a sign against a large R5 beta-cleavage share, on spent data; an item 7 side effect to register |
| MF-30 (PA in glucose + Phe, -1.7 to -3.0 dex, a lag shape) | R2 | R2's piece (1-DG + amino acid Strecker) raises early PA: **right sign**. But the same piece raises PA in M3's fed ARP-Phe argon pot, where every M3 rival already reads ~5x high (MF-29): **wrong sign** there | no net support; R2's Strecker band inherits M3's undecided route-A law (sec. 5.2) |
| MF-34 (pyrazines ~2.9 dex low) | part A (amine identity); R5 | part A's Lys class, at Yu 2017's measured Lys / Gly MGO factor (2 - 7x at fixed glucose and pH), moves the deficit by +0.3 - 0.9 dex: **right sign, at most a third of the gap**. R5 adds C5 dicarbonyls and formic acid but takes glucosone from glyoxal: sign unclear | part A's prediction is owned by M6 (CK-51). Yu 2017 is spent (sec. 2.10), so in M6 it can be a band, never decisive |

**Verdict:** no M4 rival predicts a G8 member beyond its own misfit in the right direction, at a size the cited
constants fix. The one clean prediction is part A's Lys factor on MF-34: right sign, a third of the gap at most, owned by
M6. The checks owned by M4 are CK-41 (done) and nothing else.

## 3. The rivals

| id | part G: glucose + Cys -> MFT | part E: furfural supply vs pH (deferred) | part A: amines |
|---|---|---|---|
| R0 | **current engine**, scored twice: pinned at B60's constants, and refitted on its own coordinates, on the same rows as every rival. ~0 MFT (refused); FFT through `r_glc_fur` | `acid` = [H+] | one amine, folded |
| R1 | **the claim, C-1-retained form**: Glc + Cys -> Cys-hexose Amadori -> 1-DG + Cys -> 1-DPO + C1 (C-6 lost; MFT's methyl = glucose C-1) -> the pentose steps (NF and DDP shares as M2 sets them) | amine-catalysed 1,2-enolisation, flat or weakly rising, plus an acid limb | A-R1: classes with own constants; competition declared; catalysis of dicarbonyl loss separate from consumption |
| R2 | the same supply, entering **below norfuraneol**: 1-DG + amino acid (Strecker / reduction) -> a 1,4-dideoxy C6 (Huber 1990, CF-099) -> a 1,4-dideoxy C5 (DDP) -> MFT | M1's C1: the supply carries de Roos Fig. 4's measured U-shape (flat pH 3 - 5, rising above 5.5) and the split's pH law; no separate E8 term | -- |
| R3 | **declared-failing control**: C2 + C3 recombination (B47's opened `k_glc_ha`), carrying Hofmann Table 6's measured HA + MGO -> NF and MGO's banded loss law (sec. 5.2) | the supply law is right; the steepness sits downstream (furfural loss, the HS- branch) | A-R3: one constant per alpha-amino group, whatever the carrier |
| R4 | **R4-as-written, the band-now record only**: the entry alone on the existing amine-free 1-DG supply (the red team's band-now fix), with no Cys-hexose Amadori pool. Falsified by D5 | -- | -- |
| **R4'** | **new 2026-11-06 (premortem FM3)**: the Cys-hexose Amadori pool exists, but **only as a cysteine sink**, returning cysteine and feeding no C5. The C5 entry into 1-DPO is fed by **cysteine-catalysed or amine-free enolisation** of glucose / 1-DG, not through the pool. **Decider: arm 3** (MFT per mol of fed Amadori compound consumed, against the glucose + Cys pot at matched conditions: R4' predicts the fed Amadori pot gives far less MFT per mol consumed than the sugar pot; R1 / R5 predict comparable or more) | -- | -- |
| **R5** | **the claim, C-1-lost form** (new 2026-11-06; Glomb 2010's oxidative route): Glc + Cys -> Cys-hexose Amadori (as R1) -> glucosone + Cys (B21's `r_ama_g`, oxidative) -> hydrolytic beta-cleavage: formic acid (C-1) + the C5 1,2-enediol -> 3-DPO / 1-DPO (C2 - C6; MFT's methyl = glucose C-2) -> the shared pentose steps | (its 3-DPO feeds furfural / FFT with C-1 lost, as Tressl 1994's furfuryl alcohol) | -- |

**The claim is {R1, R5}** (ruling 4). Item 6 scores R1 and R5 separately and reports **R1 vs R5 as a pair that
cannot be decided on disk** (sec. 4.3). On disk, part G can record only falsifications:
- "C2 + C3 as the hexose MFT route is falsified (R0's route and R3)";
- "R4-as-written is falsified";
- "no rival credited" (ruling 2). R0 is not displaced by any credited rival.

**Why R3 and R4 are kept as declared-failing controls.** Each is already falsified by a constant-free prediction on
spent or seen data (ruling 2: falsify, never credit). They are scored and reported in item 6 like every rival, but
cannot be favoured. R4' (sec. 3 table) is a live rival, not a control.
- **R3** fails a measured share. In both hexose pots that label MFT (Du 2023, Zhao 2019b; DEV for E1), C2 + C3 would
  leave ~50 % mixed isotopomers in a 1:1 pool, and the measured mixed share is <= ~1 %. It stays because it is the only
  rival whose MGO law Kanzler moves materially, and because the comparison should show the rejection in the record.
  - Caveat, sec. 4.5: the rejection is a transfer from 140 - 150 C, pH 6.5 pots to Hofmann's pH 3 / 5 / 7 pots.
    Without Schieberle 2003, a minor R3 share at pH 7 in Hofmann's own pots is not excluded. Item 6 reports a mixture
    R1 + f R3 with f free in [0, 0.1] as a sensitivity.
- **R4** fails structurally. It has no Cys-hexose Amadori compound, so it cannot carry the pool measured at 90 C (de
  Roos Table II) and at 145 C (Du 2024b, ~62 % of the cysteine within 10 min). That is known without a run. It stays as
  the record of the band-now fix (ENGINE_VS_LITERATURE sec. 5) only. The pool-as-sink variant is R4'.

**Notes on the rivals.**
- **R1 and R2 keep glucose C-1 as MFT's methyl; R5 keeps C-2.** A positional label on MFT separates them, unless the
  label equilibrates between C-1 and C-6 before cleavage, as Tressl 1999 Fig. 14 draws for the 3-DG side.
- **O2 does not single out R5 (FM5).** Air also drains 1-DG: x3.4 lower at 50 C (Gobert 2009), x4.4 lower in Glomb
  2010. An air / N2 contrast therefore moves R1 too, in the opposite direction.
  - The generators carry O2 on both exits.
  - The cleaner splits are arm 2's fed 1-DG vs fed glucosone under N2, and a methyl-locating label ([2-13C] or
    [1,2-13C2]-glucose; sec. 4.4).
- **The R1 / R5 mix may move with temperature (FM4).** Glomb's C-1 loss is at 50 C and Tressl's absence at 160 C, so
  the power check's generators include a temperature-dependent R1 / R5 mixture (sec. 5.4).
- **R4' and R1 / R5 share the pool and differ in what feeds 1-DPO.** D1 and D2 do not separate them (all three reach
  the same 1-DPO). Arm 3's MFT yield per mol of fed Amadori compound consumed does.
- **Corrected 2026-11-06.** The refresh said a C-1-lost route to MFT "is not proposed by any source on disk". Glomb
  2010 proposes and labels one: a C-1-lost 1-DPO, already methyl-bearing, with no reduction needed. R5 is that route.
  The decisive-measurements entry `glucose_mft_carbon_origin` still lists it as "unlisted"; its update is proposed at
  the end of this file, to be applied at the freeze by that file's one writer.
- **R5's O2.** Hofmann's aqueous pots are 100 mL in a 200 mL autoclave, with air headspace by default (no inert gas
  stated). [dry]
  - The headspace holds ~0.86 mmol O2. Glucose MFT is ~0.017 umol (1.9 ug), so MFT is not O2-limited by stoichiometry.
  - The 3.3 mmol cysteine can consume up to ~0.83 mmol O2 (4 RSH + O2 -> 2 RSSR + 2 H2O). Whether glucosone formation
    keeps up with that is a rate question no source on disk settles.
  - R5 runs at the shipped air factor (1.0) as primary, with a required-agreement sensitivity in which O2 is gone once
    the cysteine budget closes it.
- **R3 corrected (2026-11-05).** The 2026-11-04 case read R3 as "no NF from fragments". Hofmann Table 6 measures HA +
  MGO -> NF at 0.78 mol % (pH 5), more NF than furfural, so R3 predicts fragment NF whose pH law rises as in Table 6.
  The engine has no HA + MGO -> NF step today; R3 gets one at Table 6's yields, used as a fixed measured constant.
- **Part A.** A-R1 (classes) against A-R3 (one constant per alpha-amine). Luo 2022 (sec. 2.6) finds peptide /
  amino-acid ARP ratios near 1 in competition at 80 - 100 C, pH 3 - 7. That is A-R3's prediction there, against Xia
  2022's x5 - 7 in separate pots. Yu 2017 finds within-class spread (Lys / Gly x2 - 7) as large as the between-class
  contrasts. D7 is therefore not one-directional (sec. 4.2).

## 4. The decisive data

### 4.1 What cannot decide here

- **D1, D2 and D3 credit nothing (ruling 2, FM2).** D3 nominated E1. D1 was seen through R0's probe. D2's pH 3 / 7
  arms are MF-12 observations. Each may falsify a constant-free prediction and is otherwise reported as a development
  unit. The pH 3 / 7 bundles are never fitted (ruling 3).
- **Seen pots** (sec. 0) and **spent data** (sec. 2.10) can only test feasibility (development grade): Hofmann Table 1
  (D4); the ribose half of D1; Hofmann Tables 3, 4, 8, 9, 10; Zhang 2026; de Roos Fig. 4; Bolton; Meynier; Deng 2021;
  Davidek 2006a; Xia 2022; Li 2011; Steinhagen 2021.
- **D3's pots** (Du 2023, Zhao 2019b) are M2's TEST. Their route pattern is read from the published tables; no engine
  number on them before M2's item 7.
- **Band sources** (sec. 5.2) never decide, and decisive units never set a band.
- **Not simulable:** Ames 2001 (extrusion) and Cannon 2023 (no recovery, no NF). Evidence only.

### 4.2 The rows (all falsify-only: D3, D5 falsify structural predictions; D1, D2 are development units)

Directions follow from each structure. "[dry]" marks arithmetic on printed values.

| # | measurement (record) | value basis; engine status | R1 | R2 | R3 (control) | R4 (control) | R5 | separates |
|---|---|---|---|---|---|---|---|---|
| D1 | **NF / MFT in the glucose pot** (Hofmann 1998 Table 5: 7.1), with the ribose pot's 2754 as in-sample context | within-pot ratio of two printed values, no tau; **seen, falsify-only (ruling 2): a development unit, reported, never credited**; ribose cells seen (fit rows) | glucose's 1-DPO feeds NF and MFT as ribose's does, so glucose NF / MFT equals the DPO-derived ratio M2 sets. [dry] The measured NF ratio glucose / ribose is 2.5e-4 and the MFT ratio 0.096, 390x apart, so R1 needs most ribose NF to come from pentose-only routes | glucose MFT bypasses NF; NF / MFT is set only by glucose NF from other routes | NF from HA + MGO (Table 6) over MFT from HA + MP (Table 8); follows MGO's standing level, hence Kanzler's law | as R1 | **as R1** (the same 1-DPO entry and downstream shares) | R2 from R1 / R5 in development only (scale-free; nearest the intermediate). **Does not separate R1 from R5** |
| D2 | **glucose MFT pH shape** (Hofmann Table 2: 0.3 / 1.9 / 0.4 at pH 3 / 5 / 7; ribose 55.3 / 19.8 / 2.5) | within-study contrasts log10(MFT_pH / MFT_5), no tau; **seen, falsify-only (ruling 2): a development unit**; pH 3 / 7 never fitted (ruling 3) | low at pH 3 (2,3-enolisation minor in acid), falling at pH 7 as beta-cleavage takes 1-DG (Davidek 2006a) | as R1 below 5; above 5 set by the Strecker step's pH law | rising to pH 7, like fed HA + MP (Table 8) and the fragment supply (Hayashi, Yu) | set by the amine-free 1-DG supply's pH law, which no measurement fixes | set by the Amadori -> glucosone pH law, **which no measurement at >= 100 C fixes** (`k_ama_g` is licensed at pH 6.7 only); downstream as R1 | R3 from R1 / R2 by sign on the 5 -> 7 limb; the hexose shape against the pentose shape on the 3 -> 5 limb. R5 only through its declared pH sensitivity |
| D3 | **mixed isotopomers** (Du 2023, Zhao 2019b) | label share, no tau, no omega on the share (sec. 5.1); **reserved** (M2's TEST: label arithmetic only, no engine run); **spent** (nominated E1): falsify-only | intact C5 | intact C5 | ~50 % mixed in a 1:1 pool | intact C5 | intact C5 | **falsifies C2 + C3 as the hexose MFT route (R3 and R0's route), constant-free**; credits nothing; does not separate R1, R2, R4', R5 |
| D4 | **ribose / glucose MFT 10.4 and FFT 4.3** (Table 1), and the trend from > 870 at ~100 C (de Roos Table I) to 10.4 at 145 C | endpoint ratios; **seen, spent** | follows the entry's barrier. [dry] If the two labs' pots were comparable, the hexose route's apparent barrier would exceed the pentose route's by >= ~128 kJ/mol. They are not comparable (37 % water lactate vs 0.5 M phosphate) | as R1 | also follows MGO's loss law: with Ea 47 - 121 the C3 pool falls x2.5 - 33 faster at 145 C | the probe gave 45 - 78 at 145 C | as R1, plus R5's O2 budget | side evidence only |
| D5 | **the Cys-hexose Amadori pool** (de Roos Table II, 90 C; Du 2024b, 145 C; Zhang 2026 Fig. 3) with no hexose MFT at ~100 C (de Roos Table I, < 0.05 mg/mol Cys) | relative areas (II); a censored ratio > 870 (I); Du a C/C0 peak; de Roos I glucose arm **not run**, Zhang seen | the pool exists; its decay to 1-DG is the slow step at 100 C | as R1 | no Amadori step needed | no pool: fails structurally | the pool exists; its oxidative exit needs O2 | falsifies R4-as-written (structural). **R4' also carries the pool (as a cysteine sink), so D5 does not separate it** from R1 / R2 / R5; arm 3 does. Screen signatures (sec. 5.5): development grade after the screen |
| D6 | **FFT vs furfural vs pH** (Meynier Table 3; Hofmann F25#8) | within-study contrasts; seen | (part E, deferred) | | | | | deferred with part E |
| D7 | **amine class at fixed sugar**. Xia 2022: TriG-ARP x5 - 7, separate pots. Luo 2022: in-pot peptide / Gly-ARP ~1 at 80 - 100 C. de Roos Table II: Fru-Cys / Fru-Gly 1 -> 4 - 5 with dilution. Brands 2001: casein leaves glucose loss unchanged. Kanzler 2026: amines catalyse MGO loss without consumption. Yu 2017: within-class x2 - 7 | within-study ratios; Xia seen; Luo, de Roos II, Kanzler not run | (A-R1) each class its own constant and dilution law; catalysis declared | -- | (A-R3) equal per alpha-amine: Luo 2022 agrees at 80 - 100 C, pH 3 - 7; Xia and de Roos's dilution law disagree | -- | -- | **does not separate A-R1 from A-R3 on disk**: the same network's two studies point opposite ways, and within-class spread is as large as the contrast |

**The plan's decisive data, checked.**
- "Hofmann ribose / glucose ratios (MFT 10.4, FFT 4.3)": D4 here, seen and spent, side evidence.
- "The FFT pH slope": D6, deferred with part E.

### 4.3 The decisive data on disk, or their absence

- **R0: its hexose route is falsified, and no rival is credited over it (ruling 2, FM2).** R0's only hexose MFT route
  is C2 + C3 (`r_glc_c2c3`). D3's constant-free label share falsifies that route, as it falsifies R3. R0's
  no-refit residuals on D1 / D2 (-6 to -25 dex) are seen data and are reported only. **No rival is credited on disk,
  so R1 / R5 cannot ship over R0 from disk.** Recorded: "C2 + C3 as the hexose MFT route is falsified (R0's route and
  R3); no rival credited".
- **R3 (control): falsified on disk** by D3 (as above). D2's 5 -> 7 limb agrees in sign (development only).
- **R4-as-written (control): falsified on disk** by D5, structurally. **R4' is not**: it carries the pool as a
  cysteine sink (sec. 3), and its decider is arm 3.
- **The claim against R2 (the bypass): development grade only.** D1 would separate them, but D1 is seen and credits
  nothing. Under omega it
  needs >= 1.16 dex (x14) between the rivals' predictions (sec. 5.4), and its reading is reported, not credited. D2
  does not separate the two (both predict the same sign on both limbs).
- **R1 against R5 (the claim's two carbon-origin forms): cannot decide on disk, an explicit pair.**
  - D1 is identical between them by construction (one entry species).
  - D2 differs only through R5's pH law, which no measurement bands.
  - D3 is intact C5 for both.
  - The positional evidence (Tressl 1994, Glomb 2010) is labelling, the highest kind, and leans to R5's C-1 loss for
    hexose C5s in general. But it is absence without a floor, on other products, at 50 and 160 C, and Tressl's own
    laboratory disagrees with itself on the furfuryl skeleton.
- **Which glucose carbon MFT loses, a fed 1-DG, glucosone or Cys-hexose Amadori -> MFT: absent.** Sweep 2 found one
  candidate for a fed Cys-glucose ARP (P decisive ~0.1, sec. 2.9).
- **Part A:** D7 points both ways inside one network; no Maillard flavour pot varies protein-bound vs free amine at fixed
  sugar, pH and temperature. Development grade at best.

**Therefore:**
- **Part G is development grade on disk.**
  - On disk it falsifies C2 + C3 as the hexose MFT route (R0's route and R3; D3) and R4-as-written (D5).
  - It credits no rival. It reports D1 / D2 for R0, R1, R2, R4' and R5.
- **Nothing in part G ships over R0 from disk.** The decider is the designed experiment (sec. 4.4; decisive-measurements
  entry `glucose_mft_carbon_origin`):
  - arms 1 and 3 test the claim against R0 and R2;
  - arm 3 tests the claim against R4';
  - arm 2 (fed 1-DG vs fed glucosone under N2) and arm 1's labels split R1 from R5.
- **No M4 comparison runs before M2's verdict (FM8).** D1's prediction and every arm's NF / MFT read through M2's
  NF / DDP shares.
- Part A is an engine-structure requirement only (FM7, sec. 1).

### 4.4 The designed experiment (W6)

One laboratory, three arms. Throughout:
- n = 3 independent heatings per cell;
- sealed vessels with stated fill and headspace;
- pH stated on the cooled scale for every sample, and in situ where an electrode allows;
- every assay's LOD / LOQ and recovery reported;
- SIDA with labelled standards for MFT, FFT and NF.

The label, OPD and standard fixes (FM5, new 2026-11-06):
- **Quantify from unlabelled replicate pots.** In labelled pots, report isotopologue shares only. A d3-MFT internal
  standard collides with 13C-MFT's M+3 channel, so labelled pots carry no d3 standard, or use a 13C / 34S standard
  with no channel overlap.
- **OPD trapping under argon with 1 mM DTPA, with sugar + OPD and Amadori + OPD blanks** (continuous OPD makes glucosone
  from the Amadori compound; Glomb 2010). Post-incubation trapping is the primary.
- **A methyl-locating label:** [2-13C]- or [1,2-13C2]-glucose, read through an MFT fragment that keeps the methyl
  carbon. It separates C-1 (R1 / R2) from C-2 (R5) as the methyl directly, so it is not hostage to C-1 / C-6
  equilibration.

| arm | conditions | species | decides |
|---|---|---|---|
| **1. which carbon** | Hofmann's design: Cys 33 + glucose 100 mM, 0.5 M phosphate pH 3 / 5 / 7, 145 C, 20 min. Glucose as [1-13C], [6-13C], [2-13C] or [1,2-13C2], and 1:1 [13C6] / [12C6] (and unlabelled). **Each under air and under N2** (degassed, argon headspace). **Plus a second temperature, 120 C** (FM4), for [1-13C] and [2-13C] glucose with the N2 / air pair at pH 5 | MFT, NF, FFT, furfural, MP, each isotopologue-resolved, plus a methyl-locating MFT fragment | see the outcome list below |
| **2. fed C6 dicarbonyls** | synthetic 1-DG (Glomb 2000 route) **or synthetic glucosone**, 2 mM, + Na2S or Cys 6 mM, +- glycine 33 mM, 0.5 M phosphate pH 5 and 7, 145 C, 2 - 20 min, **N2 primary** (air as a secondary arm) | MFT, NF, DDP, 1-DPO and 3-DPO (OPD quinoxalines, argon + DTPA), acetic and formic acid, residual dicarbonyl | see the outcome list below; under N2 it is the cleanest R1 / R5 split (FM5) |
| **3. fed Cys-hexose Amadori** | synthetic D-fructose-L-cysteine (de Roos's preparation) 10 mM, 0.5 M phosphate pH 3 / 5 / 7, 100 and 145 C, 0 - 60 min, N2 and air; **beside a glucose + Cys pot at matched Cys, pH and temperature** | the Amadori compound, free Cys, glucosone, 1-DG and 3-DG (OPD), MFT, NF | part G's supply constants (formation, decay, Cys return, oxidative exit), the mechanism of D2's pH maximum, and **R4' vs the claim**: MFT per mol of Amadori consumed against MFT per mol of glucose-bound Cys in the sugar pot |

What each outcome decides:
- **Arm 1, [1-13C] gives a labelled MFT methyl and [6-13C] unlabelled MFT:** C-6 lost, so R1 or R2.
- **Arm 1, the reverse:** C-1 lost, so R5.
- **Arm 1, ~50 % from both:** C-1 / C-6 equilibration (Tressl 1999 Fig. 14). The methyl-locating label then decides:
  [2-13C] in the methyl means R5; [2-13C] in the ring with C-1 methyl means R1 / R2. The N2 / air contrast is
  secondary, because O2 drains 1-DG as well (FM5).
- **Arm 1 at 120 C against 145 C:** whether the R1 / R5 share moves with temperature (FM4).
- **Arm 3, the fed Amadori compound gives far less MFT per mol consumed than the glucose + Cys pot:** R4' (the pool is
  a sink). Comparable or more: R1 / R5 (the pool feeds 1-DPO).
- **Arm 1, NF's pattern equal to MFT's:** R1 or R5. Independent of it: R2. A mixed pattern: R3.
- **Arm 2, MFT forms from fed 1-DG without the amino acid, through 1-DPO and NF:** R1.
- **Arm 2, MFT needs the amino acid, with a 1,4-dideoxy quinoxaline and no NF:** R2.
- **Arm 2, MFT and C-1-lost 1-DPO form from fed glucosone under N2:** R5. **Fed 1-DG under N2 gives 1-DPO and MFT,
  and fed glucosone does not:** R1.
- **Arm 1, read again at pH 3 / 5 / 7:** D1 and D2 in a second laboratory.

Request through W6, beside the BTBU / CAAS request for MFT per glucose consumed (ENGINE_VS_LITERATURE sec. 3.1, E1).

### 4.5 What the case loses without Schieberle, Fischer & Hofmann 2003

The 2003 Weurman chapter (the first CAMOLA chapter, DFA) is still not on disk. Its role was "evidence for D3 only",
because it may share D1 / D2's laboratory and pots. Without it:
1. **R3's rejection rests on transfer.** D3's label shares come from 140 - 150 C, pH 6.5 pots in other laboratories
   (BTBU, CAAS). The chapter may hold a glucose + Cys CAMOLA on MFT in the decisive laboratory's own 145 C, 0.5 M
   phosphate design. Fragment chemistry rises with pH (Tables 6, 8; Hayashi; Yu), so a minor C2 + C3 share in
   Hofmann's pH 7 pot is not excluded. That share would bend D2's 5 -> 7 limb, the limb that separates the hexose
   shape from the rest. The R1 + f R3 sensitivity (sec. 3) is the case's answer, at the cost of one more coordinate.
2. **No DFA label on NF.** If the chapter labels NF in the glucose pot, it would show whether glucose NF is an intact C5
   (R1 / R5) or fragment-made (R3, Table 6), which is the denominator question inside D1.
3. **Nothing on which carbon.** CAMOLA ([13C6] / [12C6]) cannot tell C-1 from C-6 loss, so the chapter would not close
   the R1 / R5 gap either way.

P(it holds a glucose + Cys CAMOLA on MFT) ~0.2, unchanged. It is not a freeze blocker.

## 5. Tolerances, bands, screen and power

### 5.1 The measurement model per unit

| unit | role | value basis | reported noise | error terms in the primary | floor | density |
|---|---|---|---|---|---|---|
| D1 glucose NF / MFT (Hofmann T5) | development unit (falsify-only, ruling 2): reported, never credited | within-pot ratio | triplicates "differ by <= 10 %" (SIDA, MFT). NF method not stated: the same 10 % as a declared transfer, sensitivity x2 | omega per printed value (2 omega^2 on the ratio); no tau | MFT, NF: none printed; declared transfer 1 % of the lowest printed value in the table, required agreement over x0.1 - x10 (sec. 5.3) | exact log-normal; censored form exact |
| D1 ribose NF / MFT | context (in-sample) | within-pot ratio | as above | as above | as above | reported, not scored |
| D2 glucose MFT pH 3 / 7 against pH 5 (Hofmann T2) | development unit (falsify-only, ruling 2; never fitted, ruling 3) | two within-study contrasts sharing the pH 5 point (bivariate, correlation 0.5) | <= 10 % | omega per printed value; no tau | as D1 | exact bivariate normal in log10; censored where a rival predicts below the floor |
| D3 isotopomer shares (Du 2023, Zhao 2019b) | decisive on disk: eliminates R3 | label share | Du: no SD, two blind readers identical; Zhao: triplicate means | binomial-type share error from the printed precision (+-1 % absolute) plus a mixed-channel floor of 1 %. Omega is a level error and does not apply to a share | M+1 - M+4 channels: 1 % (the printed "< 1" cells) | exact |
| designed arms: labelled MFT (arm 1) | the decider (W6) | isotopologue share | as the request's SIDA (CV 10 %, n = 3); label read to +-1 % absolute | the **d3-MFT / M+3 collision** declared. Labelled pots report shares only; levels come from unlabelled replicate pots, and the standard's channel overlap enters as a term, or the standard is replaced | the request's LOD | exact |
| designed arms: OPD-trapped dicarbonyls (arms 2, 3) | the decider (W6) | level, ratio | as reported | an **OPD artefact term** (glucosone made from the Amadori compound under continuous trapping), measured by the sugar + OPD and Amadori + OPD blanks and subtracted with its own variance | the request's LOD | exact |
| D5 de Roos Table I (glucose < 0.05, ribose 43.5) | screen signature | censored ratio | not stated | omega per printed value on ribose; censored glucose | **0.05 mg/mol Cys (the paper's reporting limit)** | exact censored |
| D5 de Roos Table II | screen signature | relative RI areas (ratio) | response factors not given (mol ratios within ~x1.5, dossier) | omega per value plus the x1.5 response term | -- | exact |
| D5 Du 2024b pool (S5) | screen signature | peak C/C0 of the Cys-Amadori compound; post-peak loss | triplicates (figure only); the mass balance exceeds the charge by 40 % | omega per value combined with the x1.4 mass-balance term (0.444 dex in log10) | -- | exact |
| D5 Zhang 2026 MFT onset | evidence only (seen; a band source, so not screened) | ratios (RF 1) | n = 3, SD | omega per value | Zhang's printed LOD if present; else a declared transfer (1 % of the lowest printed value, required agreement) | exact censored |
| D7 class ratios (Xia, Luo 2022, de Roos, Brands, Kanzler, Yu) | part A evidence | within-study ratios; Kanzler C/C0 | as printed. Luo: no SD (0.24 - 0.34 dex per value from its scatter, dossier sec. 4.4). Brands: digitised +-3 % of span. Kanzler: +-1.5 % read, SD 0.5 - 10 %. Yu: RSD 9 % | omega per value; no tau | Kanzler: MGO <= 1 % of C0 (lowest calibrant); Yu: its printed LOQ | exact |
| R0's sub-floor predictions | all units | a predicted level against a floor is an **absolute** comparison | -- | tau_ABS 0.964 dex on the predicted level, since it is an engine level against an assay floor | the unit's floor | exact censored |
| across laboratories (D4's 100 vs 145 C trend; Bolton's bound) | side | ratios in each lab | -- | omega per value; between-lab ratio offsets are not tau (tau is absolute-only) | -- | -- |
| the MGO law (band) | band, never scored | -- | within one laboratory the MGO-alone k2 differs 0.51 dex between two nominally identical pots (Kanzler 2021 vs 2026) | the band's level and Ea ranges carry that spread (sec. 5.2) | -- | -- |

### 5.2 Bands and their sources

No band below comes from:
- Hofmann Tables 2 or 5 (seen development units, falsify-only, ruling 2);
- Du 2023 or Zhao 2019b (reserved);
- the screen's pots (de Roos Tables I and II, Hofmann Table 1, Bolton, Du 2024b).

Starts are drawn inside the bands, never set from a decisive or screened value. Some band sources have engine outputs in
other records (Li 2011, Deng 2021, Martins 2003b, Davidek 2006a, Zhang 2026; sec. 0); none is scored or screened here. A
seen source sets a centre or a fixed measured constant, and any width it does not itself print is unsourced and carries
required agreement. First-reader sources carry the rule stated at the top.

| coordinate | band | source | kind |
|---|---|---|---|
| Glc + Cys -> Cys-hexose Amadori, formation | centred so the predicted free-cysteine first-order loss at Li 2011's conditions (0.2 M each, pH 7, 100 - 130 C) matches its printed k (0.0011 /min at 100 C), x0.3 - 3; Ea 86 - 136 kJ/mol (Zhang 2026's ~101 to Li's 121, +-15) | Li 2011 (glucose + Cys itself). **FM9:** free-Cys loss also counts the thiazolidine and cystine, so it bounds net Amadori formation from **above**, not below; the band's low end therefore extends to x0.1. Du 2024b's 40 % over-closure is carried in S5's tolerance. Zhang 2026 (BTBU, units not printed: the barrier only). **One Cys-sugar Amadori class module shared with M1c** (sugar-specific factors; sec. 6) | same pot class; the width is **unsourced: required agreement** at both ends. Replaces the 2026-11-05 transfer from the trunk's Glc + Gly constant, which is kept as the required-agreement alternative |
| Cys-hexose Amadori -> 1-DG + Cys, decay | Deng 2021's Met-glucose ARP decay 0.0048 /min at 100 C (label pH 5.5), pH slope from its 7.5 arm, x0.3 - 3; **Ea 34 - 107 kJ/mol** | Deng 2021; Zhang 2026 (85.7); Martins 2005 (97 - 107); **Luo 2023 (Gly-ARP, 34 - 56; first reader)** | same class (hexose amino-acid ARP); the width unsourced: required agreement. Luo's decay is fast-then-slow; a two-pool form is not a rival here (no redesign) |
| Cys returned per Amadori decay | **0.54 - 1.0** | Martins 2003b, Davidek 2002 / 2005 (71 - 96 %); Deng 2021 (~1); Luo 2023 (0.54 - 1.0; first reader) | measured, same class |
| `r_odg_dpo` (1-DG -> 1-DPO + C1; R1) | 0 to a ceiling: the share of fed 1-DG not accounted for by the measured exits (Davidek 2006a; Voigt 2009; Smuda 2010) at their own conditions; barrier 40 - 135 kJ/mol | measured exits bound it from above only; Tressl 1994 and Glomb 2010 report no 1-DG C5 (no floor) | numerical above zero: **a verdict that needs a value off the floor does not ship without arm 2** |
| Cys-hexose Amadori -> glucosone + Cys (R5) | **0 to a ceiling (FM6)**. The ceiling is B21's shipped `k_ama_g` (Amadori -> glucosone + amine) +1.5 dex; required agreement at the floor, at the shipped value and at the ceiling. The 2026-11-06 band (the LacLys value +-0.5 dex) was too narrow for a lysine-lactulose-to-cysteine-glucose transfer. pH-flat as shipped (licensed at pH 6.7), sensitivity at the 2,3-enolisation slope (Beck / Nedvidek, ~0.44 dex/pH). O2 factor 1.0 (air), sensitivity 0 once the cysteine budget closes the O2 (sec. 3) | the shipped constant as a scale only; no Hamzalioglu value read here | numerical above zero: **unsourced, required agreement**; a verdict that needs it off the floor does not ship without arm 2's fed glucosone |
| R4''s supply (FM3) | the pool's formation, decay and Cys return as above, with **no exit to 1-DG / glucosone** (a sink that returns Cys). The 1-DPO entry is fed from glucose's amine-free enolisation (the trunk's amine-free sugar steps at their shipped values) times a cysteine-catalysis factor of 1 - 10, then `r_odg_dpo` within its band | the trunk's shipped amine-free steps; the catalysis factor has no source | unsourced: required agreement at both ends |
| O2 on the 1-DG exit (R1, and all rivals through the shared 1-DG pool) (FM5) | air drains 1-DG: x3.4 (Gobert 2009, 50 C) to x4.4 (Glomb 2010); carried as an O2-dependent 1-DG loss term, 0 to x4.4 at air, and absent under N2 | Gobert 2009; Glomb 2010 (first reader) | low-T transfer; required agreement |
| glucosone -> formic + C5 enediol -> 3-DPO / 1-DPO (R5) | 0 to a ceiling: glucosone's total measured loss at 50 C, pH 7.4 + Lys (Glomb 2010: t1/2 8 h, k ~1.4e-3 /min); barrier 40 - 135 kJ/mol; the 3-DPO : 1-DPO split as the network's pentose enediol split; on the one glucosone pool, competing with `r_g_go` | Glomb 2010 (first reader; ratio "as in glucose reactions", no numbers) | numerical above zero: **a verdict that needs a value off the floor does not ship without arm 2's fed glucosone** |
| R2's Strecker / reduction step on 1-DG | 0 to the trunk's dicarbonyl + amino-acid Strecker class constant | M3's class. M3 phase 1 ended undecidable on disk (no verdict), and M0's MF-29 finds route A's rate law, carried from fed-dicarbonyl pots, ~5x high in the argon pot, so the upper edge is suspect | class bound; required agreement |
| NF and DDP shares downstream of 1-DPO | as M2 leaves them; if M2 has no verdict, M2's rivals carried as a factor | M2 | inherited |
| **MGO loss law (`k_mgo_mel` and any MGO sink on the network)** | **Ea 47 - 121 kJ/mol**, with the level at 100 C pinned to B60's 0.0228 /min (an amine-pot level; Pfeifer's amine-free 130 C value is x7 lower); order 1 primary, order 2 a sensitivity (identified within Kanzler 2021's series); B60's Ea 20 a required-agreement sensitivity throughout | Kanzler 2026 (MGO + amino acids, 95 - 121; second reader done for CF-122); Kanzler 2021 with Bork 2022 (MGO alone, 47 - 66; first reader); Pfeifer 2013 (130 C level check; first reader) | measured class, **one laboratory lineage**; within-lab between-paper spread 0.51 dex |
| amine catalysis of MGO loss (part A) | excess k proportional to [amine]^0.75 - 0.84 (Ala); Lys x5 at 100 mM | Kanzler 2026 | one laboratory; part A evidence; no class constant fitted on flavour rows |
| R3's HA + MGO -> NF | Hofmann Table 6 yields at pH 3 / 5 / 7 | Hofmann 1998 Table 6 (not run, not scored; the same paper as D1 / D2 but other pots), a fixed measured constant | measured constant |
| R3's fragment-supply pH law (`k_glc_ha` vs pH) | **0.28 - 0.76 dex per pH unit, rising** | Hayashi 1986 (acid limb >= 0.28, C3 0.33 - 0.41); Yu 2017 (MGO standing level 0.50 - 0.76 on 5 -> 7); both first reader | measured range, standing levels not supply constants: required agreement at both ends |
| part A class constants (structure only: built at parity, off by default; FM7) | peptide / amino-acid ARP formation ratio 0.43 - 7 (Luo 2022 in-pot to Xia 2022 separate pots); Lys / Gly dicarbonyl factor 2 - 7 (Yu 2017); protein lysine from B20 | Luo 2022, Xia 2022, Yu 2017, B20 | declared bands per class; never fitted on flavour rows |
| omega | 0.419 dex; sensitivity 0.390 | H0 (EX-059) | measured |

### 5.3 Floors

- Hofmann 1998 prints no LOD for its SIDA thiols or its NF assay. The floor is a declared transfer (1 % of the lowest
  printed value in the same table), with required agreement: every D1 / D2 verdict must hold for floors from x0.1 to
  x10 of that value.
- de Roos Table I's "< 0.05 mg/mol Cys" is the paper's own reporting limit and is used as printed.
- Kanzler's lowest calibrant (1 % of C0) is its floor.
- The MFT non-detects without an LOD (Zhang & Ho 1991b, Tai 1997 / 1998, Cannon 2023 at 130 C) are evidence, not
  units.

### 5.4 Power check

**Dry arithmetic on the printed values, under omega 0.419 dex per value.** The one-degree-of-freedom threshold is 1.92
nats; a favoured rival exactly at the datum reaches it only if the other rival misses by >= 1.96 sd.
- **D1 (glucose NF / MFT, two values):** sd = 0.419 x sqrt(2) = 0.593 dex.
  - It decides between two rivals only if their predicted glucose NF / MFT differ by >= 1.16 dex (x14).
  - Suppose one rival predicted equal NF / MFT ratios in the glucose and ribose pots (the naive R1). The
    ratio-of-ratios form (four values, sd 0.838 dex) would then put the datum (-2.59 dex) 3.1 sd away, ~4.8 nats.
  - R1 does not predict that, because ribose NF has pentose-only sources.
  - **R1 vs R5: zero by construction** (identical prediction).
- **D2, 5 -> 7 limb:** measured -0.68 dex. R3's transfer from Table 8 alone is +0.07 dex, a difference of 0.74 dex =
  1.25 sd, **~0.8 nats: below the threshold**. R3's fragment-supply law (0.28 - 0.76 dex/pH) adds to that rise, while
  the MP standing level in the same pots (Table 9, seen) falls 0.35 dex from 5 to 7. So R3's net D2 prediction is not
  fixed by dry arithmetic. D3 is R3's rejection.
- **D2, 3 -> 5 limb:** measured +0.80 dex against the ribose (pentose) shape -0.45 dex. The difference is 1.25 dex =
  2.1 sd, ~2.2 nats. It rejects "the hexose route has the pentose route's pH law" marginally (one limb; the shared pH 5
  point makes the bivariate form a little weaker). That is a development reading only, under ruling 2.
- **D3:** a 1 % mixed channel against R3's ~50 %: decisive at any omega (a share, not a level).

**The power check: a reduced forced-input mini-model (SPEC; not built).** The check runs on a reduced model, not on
the full network (ruling 1 as corrected; premortem adopted 2026-11-06). It is required before the freeze. Synthetic
data only: no measured value of any unit is read, and no engine run on D3's pots (M2's TEST).
- **Precedent:** M3 phase 1's forced-input layer (EX-049 reproduced bit-exact, mass balance to ~1e-14, the
  adaptive-tempering SMC density estimator).
- **Sizing:** ~1 day of work.
- **Ordering:** it is a power check, not a comparison, so FM8 (no M4 comparison before M2's verdict) does not block it.
  M2's NF / DDP shares enter as a factor.

**Its jobs** (ruling 2: falsify, never credit):
- **(a) The on-disk falsifications hold at every band corner.**
  - R3 / R0's C2 + C3 route against D3: label arithmetic only. A C2 + C3 route in a 1:1 [13C6] / [12C6] pool leaves
    a binomial mixed share, compared with the printed <= 1 %. No engine run on Du 2023 / Zhao 2019b.
  - R4-as-written against D5: structural, no run needed.
- **(b) Development record:** D1's and D2's power, reported and never credited.
- **(c) The decider's design power.** For each arm of sec. 4.4, at n = 3 and n = 6 heatings per cell, give P(it
  separates):
  - the claim from R0;
  - the claim from R2;
  - the claim from R4';
  - R1 from R5.

  This sizes the W6 request.

**The mini-model.**
1. **Species** (one well-mixed aqueous pot, isothermal at the declared temperature):
   - glucose (forced input: its consumption clamped to a declared first-order course within Li 2011's band, or held
     at charge);
   - free Cys (state);
   - the Cys-hexose Amadori pool (state);
   - 1-DG, glucosone and 3-DG (states);
   - the C5 1,2-enediol (lumped into 3-DPO / 1-DPO at the declared split);
   - 1-DPO, NF, DDP (for R2), MFT and formic acid (states);
   - a lumped "other exits" sink per dicarbonyl (acetic acid and C4 fragments for 1-DG; glyoxal for glucosone through
     `r_g_go`'s band);
   - H2S, as a forced input proportional to Cys loss within the M1 envelope, declared.
2. **Steps per rival** (bands from sec. 5.2; no fitted constant from any pot scored here):
   - **R1:** pool -> 1-DG + Cys -> 1-DPO + C1 (C-6 lost).
   - **R5:** pool -> glucosone + Cys -> formic (C-1) + C5 enediol -> 1-DPO / 3-DPO.
   - **R2:** 1-DG + amino acid -> 1,4-dideoxy C6 -> DDP -> MFT.
   - **R4':** the pool is a Cys sink only; glucose -> 1-DG by amine-free enolisation x the Cys-catalysis factor.
   - **R0:** no entry; C2 + C3 only, at its shipped `k_glc_ha`.
   - **Shared downstream:** 1-DPO -> NF and 1-DPO -> MFT directly, NF + H2S -> MFT, all at M2's shares carried as a
     factor (each M2 rival, or M2's verdict once it exists).
   - O2 enters as a declared factor on the oxidative steps (`r_ama_g`-like) and on the 1-DG drain (FM5). It is 1 under
     air and 0 under N2, with the cysteine-budget sensitivity.
3. **Label maps.** Every carbon of glucose (C-1 ... C-6) is tracked through each step as a position vector:
   - 1-DG keeps C1 - C6 (methyl = C-1);
   - R1's 1-DPO drops C-6;
   - R5's glucosone loses C-1 to formic acid, so its 1-DPO methyl = C-2;
   - C2 + C3 pairs fragments at random from the labelled and unlabelled pools.

   MFT's carbons, and the methyl-locating fragment's, are read from the 1-DPO map. Generating variants:
   - **C-1 / C-6 equilibration before cleavage** (Tressl 1999 Fig. 14) is a variant: a 50 % swap on the hexosulose.
   - **Temperature-dependent R1 / R5 mixture (FM4):** the share of 1-DPO from R5 is f(T) = logistic in 1/T, with f
     at 145 C drawn on [0, 1] and dEa(R5 - R1) drawn on [-60, +60] kJ/mol. Arm 1's 120 C cells test whether f moves.
4. **Designs:** sec. 4.4's arms exactly (charges, buffer, pH label, temperature, time, headspace), plus Hofmann's
   aqueous arm for (b).
5. **Measurement model:** sec. 5.1 exactly:
   - omega per value (0.390 as a sensitivity);
   - SIDA CV 10 %, n = 3 or 6;
   - label shares +-1 % absolute;
   - the d3-MFT / M+3 channel collision as a declared overlap term in labelled pots (levels taken from the unlabelled
     replicate pots);
   - the OPD artefact as a blank-subtracted term with its own variance;
   - floors from the request's LODs.
6. **Draws:** Monte Carlo over every sec. 5.2 band and both ends of each required-agreement range, including:
   - the first readers' pre-dossier bands;
   - the MGO law at Ea 47 / 121 / 20, orders 1 and 2;
   - R5's `k_ama_g`-scaled band from 0 to the ceiling;
   - the O2 factors;
   - the R4' catalysis factor;
   - M2's shares.
7. **Scoring:** exact densities where closed-form (label shares, log-normal levels), else the M3 SMC estimator
   (PSIS k-hat < 0.7, ESS >= 100). Also the adequacy rule, Holm across arms, and a threshold of 1.92 nats.
8. **Checks before any result is read:**
   - mass and carbon balance per run to < 1e-10;
   - label conservation (every carbon accounted for at every step);
   - a hand-computed label table for each rival at one parameter point, reproduced exactly;
   - R0's C2 + C3 mixed share reproduced against the binomial arithmetic.
9. **Report:**
   - per arm, per pair, and per n, P(the generating rival wins by > 1.92 nats) and P("both inadequate");
   - the arms ranked by design power per heating.

   A unit or arm below 0.5 power for a pair is marked "ranks, does not decide" for that pair. **Demotions only.**
10. **Cost (~1 day):** build 3 - 4 h; checks 2 h; runs 2 - 4 h unattended (1 core, niced); write-up 1 - 2 h.
11. **Out of scope:** the full network (W3), pH in situ (the label pH is used, as each design states it), and part A.

**P(the data on disk can separate the rivals), after the rulings and the premortem.**
- C2 + C3 as the hexose MFT route falsified (R3 and R0's route) by D3: ~0.95 (the Schieberle 2003 transfer caveat).
- R4-as-written falsified by D5: ~0.95 (structural). R4' is not separable on disk (~0.02).
- The claim (R1 / R5) credited over R0: **0**. Nothing on disk can credit a rival (ruling 2).
- The claim from R2: **0** for a credited verdict. As a development reading, P(D1's power >= 0.5) ~0.45.
- R1 from R5: **~0.02** (an explicit "cannot decide on disk" pair).
- Part A: not decided here (structure only, FM7).
- P(part G ships anything from disk): **~0**.
- Expected on-disk outcome (P ~0.85): "development grade. C2 + C3 falsified as the hexose MFT route (R0's route and
  R3), R4-as-written falsified, no rival credited. Decider: the designed experiment."
- **If the designed experiment is run as amended** (all three arms, n = 3; before the mini-model sizes them):
  - P(arm 1's methyl-locating label or arm 2 under N2 separates R1 from R5) ~0.75;
  - P(arm 1 or 3 separates the claim from R0) ~0.8;
  - P(arm 3 separates the claim from R4') ~0.65;
  - P(arm 2 separates the claim from R2) ~0.6.


**Results of the reduced-model power check (2026-11-06; draft PR #101, branch `m4-reduced-power`,
`results/validation/kinetic_core_m4_reduced_power.md`, commit 6c840545; SYNTHETIC; no measured decisive, reserved or
TEST value read).** The pre-result checks passed:
- the hand-computed label table was reproduced exactly;
- every balance closed (hexose units, carbon, cysteine, sulfide and label positions all < 1e-10 mM);
- R0's 1:1 mixed share reproduced the binomial 0.5 exactly;
- step-size error was <= 0.065 dex;
- the shipped literals matched the engine within 0.5 %.

**Job (a): the on-disk falsifications, by arithmetic.**
- **C2 + C3 (R0's route and R3) against D3 is falsified.** The route leaves a mixed share of 2q(1 - q) = 0.500 in a
  1:1 pool, whatever its rate constants (0.495 - 0.500 for q = 0.45 - 0.55). D3's printed mixed cell is <= 1 %, a miss
  of **49 sd** (~1200 nats).
  - The largest C2 + C3 share of the hexose MFT still compatible is **f <= 0.059** (1.96 sd).
  - It holds at every band corner. **One assumption:** the two fragments come from independent molecules (no in-cage
    recombination of one glucose's own C2 and C3).
- **R4-as-written is falsified structurally** (D5; no Cys-hexose Amadori species). R4' is not falsified on disk.

**Job (b): D1 / D2, development only (never credited).** P(credit the generating rival), n = 3, floors x0.1 - x10:
- D1: claim vs R2 1.00; R1 vs R5 0.00; R4' vs claim 0.00.
- D2: claim vs R2 0.01 - 0.03; R1 vs R5 0.82 - 0.87; R4' vs claim 0.59 - 0.65.
- D1 + D2: claim vs R2 1.00; R1 vs R5 0.81 - 0.87; R4' vs claim 0.59 - 0.69.

These are a development record under ruling 2. They would have credited the claim over R2 had D1 / D2 been unseen.

**Job (c): the designed experiment, per arm.** Primary = P(separates), the smaller of the two directions. Required
agreement = the minimum over every declared form at both band ends, with the worst form in brackets. **DEMOTED** where
that minimum is < 0.5 ("ranks, does not decide").

| design (heatings, n = 3) | R1 vs R5 | claim vs R0 | claim vs R2 | R4' vs claim |
|---|---|---|---|---|
| arm 1 (108) | 0.99 -> 0.73 [LOD x10] | 0.96 -> 0.67 [LOD x10] | 0.97 -> 0.73 [LOD x10] | 0.06 -> 0.00 [formation x0.1] **DEMOTED** |
| arm 2 (24) | 0.99 -> 0.83 [LOD x10] | 0.97 -> 0.83 [LOD x10] | 0.98 -> 0.77 [LOD x10] | 0.00 -> 0.00 **DEMOTED** |
| arm 3 (72) | 0.86 -> 0.47 [LOD x10] **DEMOTED** | 0.88 -> 0.67 [LOD x10] | 0.16 -> 0.00 [`k_ama_g` ceiling] **DEMOTED** | 0.90 -> 0.60 [LOD x10] |
| arm 1 + arm 3 (180) | 0.99 -> 0.73 | 0.95 -> 0.83 | 0.97 -> 0.73 | 0.91 -> 0.67 |
| **arm 2 + arm 3 (96)** | **0.99 -> 0.87** | **0.99 -> 0.87** | **0.98 -> 0.77** | **0.90 -> 0.70** |
| arm 1 + arm 2 + arm 3 (204) | 0.99 -> 0.87 | 1.00 -> 0.90 | 0.99 -> 0.80 | 0.89 -> 0.67 |

At n = 6 the minimum design's worst pair falls from 0.70 to 0.60 (the entries-floor form). So **n = 6 adds nothing**:
the limit is the band, not the replicate noise. The full per-form tables, the N2 residual-O2 form, the per-truth table
and the density diagnostics are in the results file.

**THE MINIMUM DECIDING DESIGN: arm 2 + arm 3 at n = 3, 96 heatings; worst pair 0.70 under required agreement.**
- It decides all four pairs: claim vs R0, claim vs R2, R1 vs R5, R4' vs claim.
- Arm 1 (the positional labels, 108 heatings) is not needed for any pair once arms 2 and 3 run. It adds at most
  +0.03 on the weakest pair.

**Per-arm demotions:**
- arms 1 and 2 are demoted for **R4' vs claim** (only arm 3's fed-Amadori denominator separates them);
- arm 3 is demoted for **claim vs R2** (the `k_ama_g` ceiling form) and, at n = 3, for **R1 vs R5** under LOD x10
  (0.47).

**Declared limits of the check** (from the results file):
- the densities are self-normalised importance estimates that mostly fail the case's k-hat / ESS rule; M3's SMC was not
  run;
- the verdicts rest on structural |d| far beyond the Monte-Carlo error, plus a half-set verdict-flip check (arm 3's
  claim vs R2 flips 25 %: read as uncertain);
- the forms ran at R = 30 (half-width up to ~0.18 at 0.5);
- "f moves with temperature" was not tested separately;
- the MGO law was not a coordinate (R3 is job (a) only).

**A rule problem the check found: the adequacy escape credits a pure form when the truth is a mixture.** The adequacy
rule's escape ("or the rivals' means differ by more than omega") lets a unit credit one pure form when the truth is the
FM4 R1 / R5 mixture (f* ~ U[0, 1]).
- Arm 2 credits R5 at 0.97; arm 1 credits R1 at 0.94. Both are false credits.
- The "both inadequate" reading reaches 0.95 - 1.00 for every design but arm 3, so the adequacy test itself catches
  the mixture. Only the means-differ escape overrides it.
- Under **strict adequacy** (no escape) the minimum design's primary P(separates) at n = 3 falls to:

  | pair | strict adequacy, n = 3 |
  |---|---|
  | R1 vs R5 | 0.94 |
  | claim vs R0 | 0.75 |
  | claim vs R2 | 0.75 |
  | R4' vs claim | 0.76 |

- **The pending judge revision (PR #99: adequacy required on every decisive unit) removes the escape.** When it is
  adopted, this case's power numbers (job (c) and the minimum design) **must be re-checked under the revised judge**
  before the W6 request is sent. Until then, the minimum design is sized under the current rule, with the strict
  numbers above as the expected floor.
- A mixture truth must read "both pure forms inadequate; the mixture" in item 6, not a pure-form credit.

### 5.5 The mechanism-signature screen (before any fit; demotions only)

Plug-in runs at band centres and corners, on the one network, on signatures from pots that are already seen or that
this case gives up as decisive:

| signature | pot | passes if |
|---|---|---|
| S1 the hexose route shut at ~100 C | de Roos Table I (glucose < 0.05 vs ribose 43.5 mg/mol Cys; lactate pH 5, label) | predicted ribose / glucose MFT > 870 (censored) |
| S2 a Cys-hexose Amadori pool forms in situ at 90 C | de Roos Table II (Fru-Cys / Fru-Gly 0.9 - 1.1 at 15 % water) | the pool's predicted ratio within the x1.5 response term and omega |
| S3 consistency with the shipped T1 ratio | Hofmann Table 1 (ribose / glucose MFT 10.4, pH 5) | within omega on the ratio |
| S4 Bolton's bound | glucose bouillon, 120 C | predicted Cys-derived MFT <= the bound |
| **S5 the pool at 145 C** (new 2026-11-06) | Du 2024b's control (Cys + glucose 200 / 200 mM, 0.2 M phosphate pH 6.5, 145 C) | predicted peak Cys-Amadori / Cys0 between 5 and 15 min within 0.444 dex of 0.62 (i.e. 0.22 - 1.0), and a post-peak loss >= 0.06 /min |

Zhang 2026's MFT onset (absent before 90 min) and Li 2011 are **not** screen signatures: they set bands (sec. 5.2), and
one paper's pots cannot both set a band and screen it.

After the screen, a rival that fails a signature inside all its bands is demoted (reported, not dropped). None is
added, promoted or redesigned around the screen. The screen's pots are then declared seen; D1 and D2 (spent development units) are untouched.

**Stop rule.** If no rival among R1, R2, R4' and R5 passes S1, S2 and S5 inside its bands within ~15 h, part G reads "hexose
entry not identifiable on disk" and the designed experiment goes out.

## 6. What M4 needs from the one network (W3)

- **Species:** Glc, Fru, Cys, the Cys-hexose Amadori compound (new), glucosone (one pool shared with B21's trunk
  steps), 1-DG shared by the hexose and sulfur modules, 1-DPO, 3-DPO, DDP, NF, MFT, furfural, FFT, HA, MGO, MP,
  formic acid; for R2 a 1,4-dideoxy C6. Today 1-DG and glucosone live in the trunk and 1-DPO in the sulfur lane, which
  do not compose.
- **One MGO pool** with every sink declared (`r_mgo_mel` with a banded temperature law and order, `r_mgo_mp`,
  `r_cys_actz`, `r_mgo_ak`), and an HA + MGO -> NF step for R3 (Hofmann Table 6).
- **An O2 dose term** on the oxidative steps (`r_ama_g` and R5's exit), or at least the air factor with a
  cysteine-budget sensitivity. Today `ATMOSPHERE_FACTORS` = {air: 1.0}; this is the vessel / headspace operator of M0's
  G9.
- **ONE Cys-sugar Amadori class module, shared with M1c (FM9).** The cysteine-pentose (M1c's `r_cama_form`,
  `CYS_AMADORI_ROUTE`) and cysteine-hexose Amadori compounds are one class.
  - It has sugar-specific factors on formation, decay, Cys return and exits, and one pKa table.
  - Its constants are judged once, in the joint refit.
  - M4's bands still never come from M1c's fitted constants (sec. 2.11).
- **Part A is structure only (FM7):** built in W3 at parity with today's engine, and off by default. Switching it on
  is a separate part-A case.
- **Amine classes as species, not a convention:**
  - free alpha-amino acid (with identity factors), cysteine, peptide (GSH, diglycine) and protein-bound lysine (one
    key, A2b);
  - each declares which sugars it competes for, so `compose` refuses undeclared competition;
  - a separate catalytic term: amine catalysis of dicarbonyl loss without consumption (Kanzler 2026).
- **The competitor amines' Amadori chains (M0a P1b):** Gln, Lys and Ala (and every class above) condense into an
  Amadori species that decays through the same pH-lawed steps. The acrylamide lane's one-step amine + glucose ->
  melanoidin terms (`a_*_glc`) are replaced, or kept only as a declared sink beside the chain.
- **One pH module across lanes (M0a P1 - P7):** the trunk's Amadori-decay and 3-DG-exit laws act wherever their
  steps run. It also covers:
  - the envelope's band in every lane (P5);
  - the JAX mirror's pH mode (P6);
  - the final pH where a programme declares one (P7);
  - pentose + non-glycine amine pots, answered instead of refused (P3).
- **One pot with all of them:** Glc + Cys + Gly + protein in one integration (today the lanes refuse to combine).
- **In-situ pH on the hexose Amadori and enolisation steps,** with the Cys-hexose Amadori compound's pKa in the one
  titratable table; the cooled-scale observable for label-pH comparisons; and a per-step hook for M5's buffer factor
  (every Hofmann pot is 0.5 M phosphate).
- **Isotopologue pools:** labelled and unlabelled copies of the sugar carbon, position-resolved at C-1, C-2 and C-6 (or
  atom tracking), so D3 and arm 1 are observables, not hand calculations.
- **The forced-input driver** (A2 contract sec. 2): 1-DG / 3-DG courses (Deng 2021) clamped while the entry is tested.

## Open questions for the build

1. **Order with M1 and M2.** Part G feeds the pentose steps whose NF / DDP shares M2 refits and whose supply M1's C1
   replaces. **No M4 comparison runs before M2's verdict (FM8).** Only the reduced-model power check (synthetic)
   carries M2's rivals as a factor.
2. **The Hofmann glucose pH 3 / 7 arms. Resolved (ruling 3).** The older hold-out declaration binds against
   fitting. No rival's constants are fitted on them; they are scored only, as falsify-only development units (ruling 2).
3. **The Cys-hexose Amadori step's constants.** Li 2011 now centres the formation (free-Cys loss) and Du 2024b screens
   it at 145 C. The multiplier width stays unsourced (required agreement).
4. **E8 after M1.** If C1 ships, does part E survive for `r_glc_fur` and `r_arp_tdp` alone? Decide before building the
   amine-catalysed term. Schirle-Keller's 147 kJ/mol furfural law (glucose + Cys, unbuffered) bears on `r_glc_fur`;
   so does R5's C-1-lost 3-DPO, a second hexose furfural supply.
5. **Part A's constants** come from one laboratory per class (Xia, Luo, de Roos, Brands, Tang, Kanzler, Yu). They are
   declared as bands with a per-class basis; no class constant is fitted on flavour rows.
6. **Second readers.** Kanzler 2021, Pfeifer 2013, Luo 2022 / 2023, Yu 2017, Hayashi 1986, Glomb 2010, Du 2024b and
   Tressl 1999 are first-reader dossiers. Until second readers, their bands are required-agreement sensitivities (rule at
   the top).
7. **The held-dossier screen** was done by one reader on headers and keyword counts (sec. 2.9), not by two blinded
   readers. The review decides whether that is enough.
8. **W6 request:** the BTBU / CAAS MFT yields per glucose consumed would give D4 a second laboratory. The designed
   experiment (sec. 4.4) would decide R1 vs R2 vs R5.
9. **D1 and D2. Resolved (ruling 2).** They are the E1 misfit that nominated the claim, so they are spent: scored
   development units, reported and never credited. Part G has no credited unit on disk. Its verdict is development
   grade until the designed experiment, or the sweep-2 fed-ARP candidate (only if its laboratory and measurements
   qualify; sec. 2.9), supplies one.

## What should block the freeze (for the review and the coordinator)

**Superseded by the freeze (2026-11-06).** Kept as the record of what was settled before it. The coordinator's rulings
and the adopted premortem (FM1 - FM9) settled the first rounds. Then:
1. **Done:** the reduced forced-input mini-model power check ran (sec. 5.4 results, PR #101) and its demotions are
   written in. The full network (W3) is needed for the M4 build, not for the power check.
2. **Resolved:**
   - the spent-data rule (ruling 2 as corrected by FM2: spent or seen data falsify constant-free predictions, never
     credit);
   - the pH 3 / 7 bundles are never fitted (ruling 3);
   - the claim is restated with R1 vs R5 as a cannot-decide-on-disk pair (ruling 4);
   - R4' is added (FM3);
   - R5's band is widened (FM6);
   - part A is structure only (FM7);
   - no comparison before M2 (FM8);
   - the reservoir notes and the one class module (FM9).
3. **For the review.** Test the R5 form:
   - whether glucosone from the Cys-hexose Amadori compound at 145 C is plausible while cysteine scavenges the O2;
   - whether the methyl-locating fragment of MFT is measurable (a fragment keeping the 2-methyl carbon) in the
     requesting laboratory's MS.
4. **Not blockers, but the owner should know:**
   - On disk, part G credits no rival and ships nothing over R0. It records only falsifications: C2 + C3 as the hexose
     MFT route (R0's route and R3), and R4-as-written.
   - The decider is the designed experiment, which comes last by the owner's ruling (2026-11-05).
   - The M4 build is still worth its network pieces (M5 / M6 reuse them, as does the one Cys-sugar Amadori module
     shared with M1c), and the mini-model sizes the W6 request.

## Proposed update to `docs/reference/decisive_measurements.yml` (entry `glucose_mft_carbon_origin`)

To be applied at the freeze by that file's one writer; not applied here. The same text is in PR #90's body.

**Rivals.**
- Replace "unlisted: a route through a C-1-lost C5 ..." with "R5: the oxidative route, Cys-hexose Amadori ->
  glucosone -> beta-cleavage (C-1 lost as formic acid) -> 1-DPO -> the pentose steps (Glomb 2010)".
- Restate R1 and R5 as the claim's two carbon-origin forms ("a hexose C5 entry into 1-DPO; C-1 retained vs lost"), an
  explicit pair that cannot be decided on disk.
- Add "R4': the Cys-hexose pool exists as a cysteine sink only; the C5 entry is fed by cysteine-catalysed or amine-free
  enolisation, not through the pool". Its decider is arm 3.
- Keep R4-as-written as the band-now record, falsified by D5.
- Record "C2 + C3 as the hexose MFT route is falsified (R0's route and R3) by D3's label share; no rival credited".

**`why_not_on_disk`.**
- The spent-data rule: spent or seen data may falsify a constant-free prediction, never credit a rival.
  - D3 (BTBU / CAAS label shares) nominated E1.
  - D1 (Hofmann Table 5) was seen through R0's probe.
  - D2's pH 3 / 7 arms are MF-12 observations.
- Glomb 2010 (Halle, 50 C) labels the hexose 1-DPO as C-1-lost via glucosone; nothing at >= 100 C measures it.
- O2 drains 1-DG too (Gobert 2009, Glomb 2010), so air / N2 does not single out R5.

**Arm 1.**
- Add [2-13C] or [1,2-13C2]-glucose, read through a methyl-locating MFT fragment, beside [1-13C], [6-13C] and the 1:1
  [13C6] / [12C6].
- Air / N2 sub-arms (o2: "air and N2 arms, degassed with argon headspace").
- **A second temperature, 120 C**, for the [1-13C] and [2-13C] glucose cells with the N2 / air pair at pH 5
  (temperatures_C: [120, 145]).

**Arm 2.**
- Synthetic glucosone beside 1-DG, N2 primary, with 3-DPO and 1-DPO measured.
- "Fed 1-DG vs fed glucosone under N2" is the cleanest R1 / R5 split.

**Arm 3.**
- Add a glucose + Cys pot at matched Cys, pH and temperature.
- The outcome "MFT per mol of fed Amadori consumed, far below the sugar pot's" means R4'; comparable or more means
  R1 / R5.

**Measurement requirements (all arms).**
- Levels come from unlabelled replicate pots. Labelled pots report shares only, with no d3-MFT standard: d3 collides
  with 13C-MFT's M+3 channel, so use a non-overlapping standard.
- OPD trapping under argon with 1 mM DTPA, post-incubation primary, with sugar + OPD and Amadori + OPD blanks (the
  continuous-trapping glucosone artefact).
- LOD / LOQ and recovery per assay.

**Sizing (from the reduced-model power check, PR #101; added at the freeze).**
- **Minimum deciding design: arm 2 + arm 3, n = 3, 96 heatings.**
  - It decides claim vs R0, claim vs R2, R1 vs R5 and R4' vs claim; worst pair 0.70 under required agreement.
  - n = 6 adds nothing.
  - Arm 1 (108 heatings) is optional: it adds <= 0.03 on the weakest pair.
- **Per-arm demotions** ("ranks, does not decide" for that pair):
  - arms 1 and 2 for R4' vs claim;
  - arm 3 for claim vs R2;
  - arm 3 for R1 vs R5 at n = 3 when the LOD is 10x worse than declared.
- **Caveat:** sized under the current adequacy rule. Under strict adequacy (the PR #99 judge revision) the design's
  primary power at n = 3 falls to 0.75 - 0.94. Re-check before the request is sent.
- **Set `measurement.arms`:** arm 2 and arm 3 required, arm 1 optional; replicates n = 3.

**Outcomes.**
- "C-1 lost" (unlabelled methyl from [1-13C], labelled from [2-13C]) means R5.
- "C-1 retained" means R1 / R2.
- "~50 % from both positional labels" means C-1 / C-6 equilibration, and the methyl-locating label decides.
- A shift between 120 and 145 C means a temperature-dependent R1 / R5 mixture.

## Revision log

- **2026-11-04 (v1).** Items 1 - 5 written before any result; no engine run.
- **2026-11-05 (refresh, branch `m4-case-refresh`, 8660c005), not frozen.** Still before any result; no engine run on
  any pot; the sweep read titles and abstracts only.
  - This week's rules applied: value basis per unit, tau on absolute levels only, omega 0.419, floors from assays,
    exact densities, the adequacy rule, R0 pinned at B60, bands and starts from unscored sources with required
    agreement for unsourced widths, the screen before any fit (demotions only), the power check, pH scales named.
  - Sec. 0 declares which pots have been seen under the engine. Hofmann Table 9 and Zhang 2026 moved to seen. D1's
    and D2's glucose pots are declared run under R0 only.
  - Kanzler 2026 added (fed MGO; the MGO sink's Ea vs `k_mgo_mel`'s 20 at its floor).
  - Held sources added: Tressl 1994, Voigt 2009 / 2010, Smuda 2010, Du 2025, Zhang & Ho 1991b, Martins 2010,
    Schirle-Keller 1992. Hofmann 1997, Madruga 1995 and Liu 2021 left out with reasons.
  - Hofmann Table 8 corrected (HA + MP); Table 6 (HA + MGO -> NF) added; R3's D1 prediction corrected.
  - Citation sweep 1 (12 seeds, 634 works, 511 unseen, 13 candidates with roles fixed before reading).
  - D1's primary is the glucose pot's NF / MFT; the designed experiment named.

### Changelog: 2026-11-06, READY TO FREEZE (items 1 - 5), not frozen

Still before any result of the change; no engine run on any pot. Started from the refresh (8660c005), the newer text.
- **The 14 sweep dossiers folded in by kind** (PR #76; all decisive roles NOT MET; sec. 2.9 table). New rows: Tressl
  1999 and Glomb 2010 (labelling); Schieberle 2005 (labelling, by analogy); Kanzler 2021, Pfeifer 2013, Luo 2023 and
  Glomb's half-lives (fed); Glomb Table I (trapping); Du 2024b and Cannon 2023 (intermediates); Yu 2017 (mass balance);
  Yeo 1991, Hayashi 1986, Yu 2017 and Luo 2022 (within-study). El-Ghorab 2010 and Fanun 2001 left out with reasons.
- **Tressl re-checked against both Tressl records.** The case's reading of Tressl 1994 was correct for what Tressl
  1994 prints: furfuryl alcohol, 2-hydroxymethylthiophene and the thiolanones are C2 - C6. The row now says so and
  adds that the C6 furfurylthiols keep C-1. Tressl 1999 Fig. 14 (same laboratory, no conditions, no citation) puts the
  same furfuryl skeleton at 50 % C-1-labelled. Neither share is used for that skeleton. Arm 1's interpretation now
  allows for C-1 / C-6 equilibration.
- **Glomb 2010 contradicted a line in the case** ("a route through a C-1-lost C5 ... is not proposed by any source on
  disk"). Corrected. The route is now rival **R5** (oxidative glucosone, C-1 lost, MFT's methyl from C-2), with bands
  (sec. 5.2), network needs (O2 term, glucosone pool, C-2 isotopologue; sec. 6) and two design changes (N2 / air on
  arm 1; fed glucosone on arm 2). D1 cannot separate R1 from R5; no unit on disk can.
- **Kanzler 2021 vs 2026 (x3.2 within-lab k2).** The MGO-alone Ea basis falls from ~95 - 113 to 47 - 66 kJ/mol. The
  MGO band is now Ea 47 - 121, with order 2 identified as a sensitivity and Ea 20 kept as required agreement. Pfeifer
  2013 is the same lineage, not a second laboratory. CF-122 already carries both notes.
- **Schieberle 2005's HDMF Table 3** repeats Hofmann 1998 cells. It is not used; only the primary (Hofmann Table 6)
  enters. Schieberle 2003 is still not on disk; sec. 4.5 states what that costs.
- **Rivals R3 and R4 are now declared-failing controls** (R3: D3's measured share; R4: no pool, against de Roos Table
  II and Du 2024b). They are scored and reported, but cannot be favoured.
- **Part A:** Luo 2022 points the other way from Xia 2022 inside one network, and Yu 2017's within-class spread matches
  the between-class contrasts. D7 no longer separates A-R1 from A-R3; part A stays as declared bands.
- **Bands revised from the new sources** (sec. 5.2):
  - formation centred on Li 2011's glucose + Cys free-Cys loss (held, newly cited);
  - decay Ea widened to 34 - 107 and Cys return to 0.54 - 1.0 (Luo 2023);
  - R3's fragment pH law set to 0.28 - 0.76 dex/pH (Hayashi, Yu);
  - R2's Strecker band flagged by M3's undecided route A and MF-29.
- **Screen:** S5 (Du 2024b's pool at 145 C) added.
- **The held-dossier screen** (open question 7) was done by one reader: 199 files, 135 uncited. Ames 2001 (a second
  laboratory's glucose + Cys MFT pH series, extrusion, rising with pH) entered as evidence against the 5 -> 7 fall;
  Li 2011, Tai 1997 / 1998 and Wang 2020 entered too; others were left out with reasons.
- **Sweep 2** (7 seeds, 313 unseen) was keyed to the amended criterion ((f) - (h)). It found 5 candidates with roles
  fixed before reading; DOIs in sec. 2.9.
- **New sections:**
  - 2.10, spent data (sec. 7 item 2's rule);
  - 2.11, M0's groups G4 and G8: no M4 rival predicts a member beyond its own misfit at a size cited constants fix; part
    A's Lys factor gives MF-34 the right sign for at most a third of the gap (owned by M6);
  - 4.5, what is lost without Schieberle 2003;
  - the synthetic power check's specification (sec. 5.4);
  - P(the data on disk separate the rivals);
  - what should block the freeze;
  - the proposed update to the decisive-measurements entry.
- The private download queue is no longer named by path.

### Changelog: 2026-11-06 (later), coordinator rulings and M0a folded in; still READY TO FREEZE, not frozen

- **Ruling 1:** the synthetic power check waits for the one JAX network (W3). Its specification (sec. 5.4) now also
  sizes the designed experiment's arms, since the decider is W6.
- **Ruling 2:** D1 and D2 (Hofmann 1998's glucose pots) are spent. They are the E1 misfit that nominated the claim.
  - They are scored development units, reported and never credited (secs. 0, 2.10, 4.1 - 4.3, 5.1).
  - Part G is development grade on disk. It can eliminate R3 (D3) and R4 (D5).
  - The claim cannot ship over R0 from disk. The decider is the designed experiment (entry `glucose_mft_carbon_origin`).
  - The odds are rewritten: P(part G ships anything from disk) ~0.
- **Ruling 3:** the older hold-out declaration on the glucose pH 3 / 7 bundles binds against fitting. They are never
  fitted; they are scored only.
- **Ruling 4:** item 1 is restated as "a hexose C5 entry into 1-DPO; carbon origin (C-1 retained vs lost) undecided".
  R1 (C-1 retained, through 1-DG) and R5 (C-1 lost, through glucosone) are its two forms, an explicit "cannot decide on
  disk" pair. The rivals table, odds, power check and proposed decisive-measurements update follow.
- **M0a (draft PR #91) is folded into part A.** It adds the bypass list P1 - P7 (the trunk pH law runs only in the
  trunk lane) and P1b (the acrylamide lane's Gln / Lys / Ala react through one amine + glucose -> melanoidin step with
  no Amadori chain).
  - Part A now requires amine classes as species with their Amadori chains and declared competition, plus one pH module
    across lanes (secs. 1, 6).
  - Part A inherits M0a's spent data.
- **Yu 2017 is listed as SEEN / spent** (sec. 2.10). It was used for R3's fragment pH band and part A's within-class
  spread, so it is credited as decisive nowhere.

### Changelog: 2026-11-06 (premortem adopted, FM1 - FM9); READY TO FREEZE as development grade, decider W6, pending the reduced-model power check

- **Status** changed to "READY TO FREEZE as development grade, decider W6: pending the reduced-model power check".
- **FM2 corrects ruling 2's basis.** D3 (the BTBU / CAAS label shares) is what nominated E1. D1 was only seen through
  R0's probe, and D2's pH 3 / 7 arms are MF-12 observations.
  - **The rule:** spent or seen data may falsify a rival whose prediction is constant-free, never credit one. It
    applies to R0's route as to R3.
  - **Recorded:** "C2 + C3 as the hexose MFT route is falsified (R0's route and R3); no rival credited".
    R4-as-written is falsified by D5.
- **FM3:** new rival R4' (the pool as a cysteine sink; the C5 entry from cysteine-catalysed or amine-free enolisation).
  Its decider is arm 3 (MFT per mol of fed Amadori consumed against a matched glucose + Cys pot). R4-as-written is kept
  as the band-now record only.
- **FM4:** a temperature-dependent R1 / R5 mixture in the generators; arm 1 gains a 120 C cell ([1-13C], [2-13C]
  glucose, N2 / air, pH 5).
- **FM5:**
  - O2 drains 1-DG too (Gobert, Glomb), so air / N2 no longer singles out R5;
  - a methyl-locating [2-13C] / [1,2-13C2] label added;
  - arm 2 under N2 is the primary R1 / R5 split;
  - OPD under argon + DTPA with sugar / Amadori + OPD blanks;
  - levels taken from unlabelled replicate pots;
  - the OPD artefact and the d3-MFT / M+3 collision added to the measurement model.
- **FM6:** R5's Amadori -> glucosone step is banded "0 to a ceiling" (shipped `k_ama_g` + 1.5 dex), not the LacLys
  value +-0.5.
- **FM7:**
  - part A is an engine-structure requirement only (W3, at parity, off by default; switching on is a part-A case);
  - M0a's sweep candidates are carried into part A with M0a's roles (Wang Y. 2019, 10.1016/j.foodres.2019.05.001,
    among them);
  - the pH scale is resolved: in situ where a lane has a pH state, the declared label pH elsewhere;
  - Hollnagel & Kroh 1998 is reconciled: out of part G, an evidence candidate for part A.
- **FM8:** no M4 comparison runs before M2's verdict.
- **FM9:**
  - Li 2011's free-Cys loss includes thiazolidine and cystine, so it bounds Amadori formation from above (the band's
    low end goes to x0.1);
  - Du 2024b's 40 % over-closure is kept in S5;
  - one Cys-sugar Amadori class module with sugar-specific factors is shared with M1c.
- **Sec. 5.4** now specifies the reduced forced-input mini-model power check (pool -> 1-DG / glucosone -> 1-DPO ->
  NF -> MFT with carbon label maps, M3-layer style, ~1 day). It is not built. The D3-design engine runs are dropped
  (M2's TEST; label arithmetic suffices). The odds are rewritten.
- **The decisive-measurements proposal** is updated: R4', arm 3's matched pot, the label / OPD / standard fixes, and the
  second temperature. It is in the PR body too.

### Changelog: 2026-11-06 (later still), two candidates obtained and read; status unchanged

- The roles of Lei 2024 (10.1016/j.foodres.2024.114454) and Wang Y. 2019 (10.1016/j.foodres.2019.05.001) were
  confirmed and committed (a29535db) before either paper or its SI was opened.
- Dossiers `lei2024_extraction.md` and `wang2019_extraction.md` were added, and both SI files read under the same
  roles.
- Lei 2024: the decisive role is **NOT MET** (no NF, no recovery, dry; SI is DFT only, excluded). It is weak evidence
  for D5.
- Wang 2019: part A's role is **MET IN FORM, CONFOUNDED** (pH x 0.4 M acetate; the cooled pH converges). It goes to the
  later part-A case.
- Neither changes "freeze as development grade, decider W6" or "build only the W3 structure plus the reduced-model
  power check".

### Changelog: 2026-11-06, FROZEN (items 1 - 5)

- Sec. 5.4 records the reduced-model power check (PR #101, commit 6c840545).
  - **Job (a):** C2 + C3 is falsified against D3 by arithmetic (mixed share 0.5 vs <= 1 %, 49 sd). The C2 + C3 share
    of hexose MFT is f <= 0.059. Assumption: the fragments come from independent molecules. R4-as-written is falsified
    structurally.
  - **Job (b):** D1 / D2 are development only.
  - **Job (c):** the per-arm table, with demotions:
    - arms 1 / 2 for R4' vs claim;
    - arm 3 for claim vs R2;
    - arm 3 for R1 vs R5 at n = 3 under LOD x10.
  - **Minimum deciding design:** arm 2 + arm 3, n = 3, 96 heatings (worst pair 0.70 under required agreement);
    n = 6 adds nothing.
- **The rule problem is recorded.** The adequacy rule's "means differ by more than omega" escape credits a pure form
  when the truth is the FM4 mixture. The pending judge revision (PR #99: adequacy required on every decisive unit)
  removes it. The case's power numbers must be re-checked under the revised judge once it is adopted.
- **Status FROZEN:** "development grade on disk; decider W6 (`glucose_mft_carbon_origin`), sized: arm 2 + arm 3,
  n = 3, 96 heatings". The freeze commit's sha is recorded by the commit after it.
- Added at the freeze (coordinator): item 5 tolerances follow sec. 7 item 5 as revised (branch `judge-revision`,
  455d8396; adoption pending). Before any scoring, the measurement model and the reduced power check are re-checked
  under the adopted text, which also removes the means-differ adequacy escape. No decisive data is read before then.
