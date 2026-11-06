# Case M4: glucose to MFT; protein and other amines as classes

**READY TO FREEZE (items 1 - 5), not frozen** (2026-11-06). Items 1 - 5 are complete under this week's rules. They
freeze when three things have happened, in this order:
1. the independent review (one round plus one confirmation) has answered the questions in "What should block the
   freeze" (end of this file);
2. the synthetic power check of sec. 5.4 has run on the first build of the one network (W3);
3. the demotions it orders are written in.

No engine was run for this case, on any pot. The changelog is at the end.

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
  ribose and xylose T1 rows). See open question 2.

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
- "Run under R0 only" = the incumbent's no-refit probe has a number there and no other structure has run. Such a pot
  stays decisive between the new rivals, because R0's numbers are refusal-scale and predate every rival's design. The
  review is asked to check this (open question 9).

| pot (record) | status | where |
|---|---|---|
| Hofmann 1998 Table 1 (= Table 2's pH 5 arms) | **seen** | ribose, xylose: B60 fit rows; fructose, glucose: `thiol_submodel.yml` rows and B8's level rows; T1 design in probe E and the red-team probe |
| Hofmann 1998 Table 2, ribose pH 3 / 7 | **seen** | `thiol_submodel.yml` rows; probes B and D in shape (30 mM design) |
| Hofmann 1998 Table 2, **glucose pH 3 / 7** | **run under R0 only** (refused in the core panel; R0's residuals in the no-refit thiol sub-model probe, summary seen) | mp hold-out bundles; `thiol_submodel.yml` rows `hofmann1998_cys_glucose_ph{3,7}` |
| Hofmann 1998 Table 5, ribose NF and furfural | **seen** | B60 fit rows |
| Hofmann 1998 Table 5, **glucose NF / furfural** | **run under R0 only** (the same probe, summary seen); rhamnose not run | `hofmann1998_cys_glucose_ph5_fa_nf` |
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
| G. hexose entry to MFT | new: Glc + Cys -> cysteine-hexose Amadori compound (D-fructose-L-cysteine) -> 1-DG + Cys released; then `r_odg_dpo` (1-DG -> 1-DPO + C1, [sulfur.py](../../../src/kinetic_core/sulfur.py) ~927, "the C-6 loss") into the shared pentose MFT steps | see below | a measured-class Cys-hexose Amadori formation and decay; `r_odg_dpo` on with a declared barrier |
| E. the furfural supply's pH law (E8) | `r_pent_tdp`, `r_arp_tdp`, `r_glc_fur` (`acid` = [H+], a decade per unit; sulfur.py ~1160 - 1175) | as stated | an amine-catalysed 1,2-enolisation term, pH-flat or weakly rising, beside a weaker acid limb. **Deferred** (plan M4 row) |
| A. amine classes | the trunk's `r_schiff` (Glc + Gly), the sulfur lane's Cys-catalysed `r_pent_dpo` / `r_pent_tdp`, B20's `r_glc_lysp` (protein lysine, trunk only) | see below | classes: free alpha-amino acid, cysteine, peptide, protein-bound lysine; each with its own condensation constant; their competition for each sugar declared; **amine catalysis of dicarbonyl loss declared separately from amine consumption** (Kanzler 2026, sec. 2.2) |

**Today, part G.** Glucose reaches MFT only through `r_glc_c2c3` (log10 k -8.59 at 145 C; Glc -> HA + MGO), then
`r_mgo_mp` (MGO + H2S -> 1-mercapto-2-propanone, MP) and `r_ha_mp_mft`. FFT forms through `r_glc_fur` (`acid`
factor). `HEXOSE_C5_ENTRY.k_per_min_100c = 0` ([parameters_sulfur.py](../../../src/kinetic_core/parameters_sulfur.py)
~623 - 641). Cysteine is not an amine on the trunk, so a Glc + Cys pot has no amine route to 1-DG.
`engine.unidentified_routes` refuses hexose-only thiol requests (engine.py ~1024 - 1043).

**Today, part A.** There is one amine: MET and PRO are folded into `Gly` (engine.py ~2560 - 2570). `LYSP` exists only
on the trunk lane. There is no peptide.

**The claim's weak point (new 2026-11-06).** The claim's specific C5 step is 1-DG -> 1-DPO + C1, which keeps glucose
C-1 as MFT's methyl. Two laboratories' positional labels speak against it, weakly:
- **Glomb 2010** (Halle; 50 C, pH 7.4, Lys, air vs de-aerated) finds the hexose C5 dicarbonyls, 1-DPO included,
  "strictly" glucose C2 - C6. They come from glucosone by hydrolytic beta-cleavage and rise with O2; none is assigned
  to 1-DG.
- **Tressl 1994** (TU Berlin; 160 C, pH 5, Cys) reports no C5 product from 1-DG; every 1-DG product keeps six carbons.

Both are absences without a floor, neither measures MFT, and neither is at Hofmann's 145 C, pH 3 - 7. So they do not
refute the claim, but they name a rival with the same entry species and the opposite label (R5, sec. 3).

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
| Hofmann 1998 Table 5 (`hofmann1998_t5_fa_nf_pairs`; dossier sec. 2.5) | Cys 33 + sugar 100 mM, 0.5 M phosphate pH 5, 145 C, 100 mL in a 200 mL autoclave: NF / MFT 2754 (ribose), 7.1 (glucose), 23.9 (rhamnose); furfural / FFT 5.58 / 2.46 / 7.75. NF and furfural method not stated (MFT, FFT by SIDA) | whether glucose MFT passes through the intermediate that ribose MFT does (D1) | DEV; ribose NF / furfural are fit rows; the glucose cells run under R0 only (sec. 0) |
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
| Hofmann 1998 Table 2 (F23#18, F25#8; dossier sec. 2.2) | 0.5 M phosphate pH 3 / 5 / 7: glucose MFT 0.3 / 1.9 / 0.4, FFT 0.7 / 2.8 / 0.6; ribose MFT 55.3 / 19.8 / 2.5, FFT 22.9 / 12.1 / 1.2 | the hexose route's pH shape against the pentose one (D2) | DEV; glucose pH 3 / 7 are Maillard-path hold-out bundles, **run under R0 only** |
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

**Held dossiers (2026-11-06 screen).** Open question 7 asked for a screen of every held dossier that names glucose
with a thiol or a deoxyosone. That screen found 199 such files; 135 were not cited in the 2026-11-05 text. Each was
screened on its header and on keyword counts by one reader (not two blinded readers, as M1c did). Those that bear on
the criterion now sit in sec. 2 or 2.8 with their reason: Ames 2001, Li 2011, Tai 1997 / 1998, Wang 2020, US 9,005,689,
Haleva-Toledo 1999, Zyzak 1995 and Huyghues-Despointes 1996. The rest are pentose-only pots, thiamine pots, meat and
other food matrices, patents on heme systems, trunk-only dicarbonyl kinetics (Goncuoglu, Gokmen, Yu 2018 / 2020,
Zhang 2020), downstream fed furanones (Zheng 1997, US 3,904,655, Munch 1997 / 1998), and reviews. Ames 2001 is the one
the case had missed that speaks against a rival (sec. 2.6).

### 2.10 Spent data (sec. 7 item 2)

**Spent: these data nominated a rival, and are development evidence for it.**

| data | spent by | spent for | consequence |
|---|---|---|---|
| Hofmann 1998 Table 1 (all arms) and the T1 design | atlas E1; probes C and E; the red-team probe; M0 CK-41 (MF-12's size); B8 - B51.1 level rows | R1, R4 (the claim and the band-now fix were nominated on it) | D4 is side evidence only; S3 is a screen signature |
| Bolton 1993 / 1994 | core panel; M0 MF-15 | -- (no M4 rival was nominated on it) | S4 only |
| Glomb 2010 Table I and its label statements | the 2026-11-05 dossier, which nominated R5 | R5 | never decisive (it is not a thiol pot anyway) |
| Steinhagen 2021 (glucose alone, 121 C) | M0 MF-23 (read for sec. 2.11) | R5's oxidative piece | sign evidence only |
| hofmann2000 Fig. 2, Tables 3 - 4 (glucose + Phe PA) | M0 CK-45 (planned), MF-30 | R2's 1-DG + amino-acid Strecker piece | sign evidence only (sec. 2.11) |
| Kanzler 2021 / 2026, Pfeifer 2013, Bork 2022 | CF-122, M0 MF-24, CK-46 / CK-47 (planned) | the MGO loss law (a band, not a rival) | band source only, never scored |
| Hofmann 1998 Tables 6, 8 | the 2026-11-05 R3 correction | R3 | fixed measured constants for R3, never scored |
| Du 2023 / Zhao 2019b route pattern (published tables, no engine run) | the 2026-11-04 case | R3's rejection | D3 rejects R3 as a measured share; no engine number on them before M2's item 7 |

**Seen but not spent.** These have engine numbers on record but nominated no rival here.
- Hofmann 1998 Table 2 glucose pH 3 / 7 and Table 5 glucose NF / MFT: run under R0 only (sec. 0). The claim was
  written on 2026-11-04 from Table 1 and the labelling; R0's probe summary was read on 2026-11-05. These stay D1 / D2
  unless the review rules that any engine number spends a pot (open question 9).
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
| MF-34 (pyrazines ~2.9 dex low) | part A (amine identity); R5 | part A's Lys class, at Yu 2017's measured Lys / Gly MGO factor (2 - 7x at fixed glucose and pH), moves the deficit by +0.3 - 0.9 dex: **right sign, at most a third of the gap**. R5 adds C5 dicarbonyls and formic acid but takes glucosone from glyoxal: sign unclear | part A's prediction is owned by M6 (CK-51). If M6 uses Yu 2017 for a Lys factor, Yu is spent for part A there |

**Verdict:** no M4 rival predicts a G8 member beyond its own misfit in the right direction, at a size the cited
constants fix. The one clean prediction is part A's Lys factor on MF-34: right sign, a third of the gap at most, owned by
M6. The checks owned by M4 are CK-41 (done) and nothing else.

## 3. The rivals

| id | part G: glucose + Cys -> MFT | part E: furfural supply vs pH (deferred) | part A: amines |
|---|---|---|---|
| R0 | **current engine**, scored twice: pinned at B60's constants, and refitted on its own coordinates, on the same rows as every rival. ~0 MFT (refused); FFT through `r_glc_fur` | `acid` = [H+] | one amine, folded |
| R1 | **the claim**: Glc + Cys -> Cys-hexose Amadori -> 1-DG + Cys -> 1-DPO + C1 (C-6 lost; MFT's methyl = glucose C-1) -> the pentose steps (NF and DDP shares as M2 sets them) | amine-catalysed 1,2-enolisation, flat or weakly rising, plus an acid limb | A-R1: classes with own constants; competition declared; catalysis of dicarbonyl loss separate from consumption |
| R2 | the same supply, entering **below norfuraneol**: 1-DG + amino acid (Strecker / reduction) -> a 1,4-dideoxy C6 (Huber 1990, CF-099) -> a 1,4-dideoxy C5 (DDP) -> MFT | M1's C1: the supply carries de Roos Fig. 4's measured U-shape (flat pH 3 - 5, rising above 5.5) and the split's pH law; no separate E8 term | -- |
| R3 | **declared-failing control**: C2 + C3 recombination (B47's opened `k_glc_ha`), carrying Hofmann Table 6's measured HA + MGO -> NF and MGO's banded loss law (sec. 5.2) | the supply law is right; the steepness sits downstream (furfural loss, the HS- branch) | A-R3: one constant per alpha-amino group, whatever the carrier |
| R4 | **declared-failing control**: the entry alone on the existing amine-free 1-DG supply (the red team's band-now fix) | -- | -- |
| **R5** | **new 2026-11-06, the oxidative C-1-loss route** (Glomb 2010): Glc + Cys -> Cys-hexose Amadori (as R1) -> glucosone + Cys (B21's `r_ama_g`, oxidative) -> hydrolytic beta-cleavage: formic acid (C-1) + the C5 1,2-enediol -> 3-DPO / 1-DPO (C2 - C6; MFT's methyl = glucose C-2) -> the shared pentose steps | (its 3-DPO feeds furfural / FFT with C-1 lost, as Tressl 1994's furfuryl alcohol) | -- |

**Why R3 and R4 are kept as declared-failing controls.** Each is already known to fail on development data. They are
scored and reported in item 6 like every rival, but cannot be favoured.
- **R3** fails a measured share. In both hexose pots that label MFT (Du 2023, Zhao 2019b; DEV for E1), C2 + C3 would
  leave ~50 % mixed isotopomers in a 1:1 pool, and the measured mixed share is <= ~1 %. It stays because it is the only
  rival whose MGO law Kanzler moves materially, and because the comparison should show the rejection in the record.
  - Caveat, sec. 4.5: the rejection is a transfer from 140 - 150 C, pH 6.5 pots to Hofmann's pH 3 / 5 / 7 pots.
    Without Schieberle 2003, a minor R3 share at pH 7 in Hofmann's own pots is not excluded. Item 6 reports a mixture
    R1 + f R3 with f free in [0, 0.1] as a sensitivity.
- **R4** fails structurally. It has no Cys-hexose Amadori compound, so it cannot carry the pool measured at 90 C (de
  Roos Table II) and at 145 C (Du 2024b, ~62 % of the cysteine within 10 min). That is known without a run. It stays as
  the record of the band-now fix (ENGINE_VS_LITERATURE sec. 5) and as the control that shows whether the amine-borne
  entry matters on D1 / D2.

**Notes on the rivals.**
- **R1 and R2 keep glucose C-1 as MFT's methyl; R5 keeps C-2.** A positional label on MFT separates them, unless the
  label equilibrates between C-1 and C-6 before cleavage, as Tressl 1999 Fig. 14 draws for the 3-DG side. Sec. 4.4,
  arm 1, adds an N2 / air contrast for that case.
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

- **Seen pots** (sec. 0) and **spent data** (sec. 2.10) can only test feasibility (development grade): Hofmann Table 1
  (D4); the ribose half of D1; Hofmann Tables 3, 4, 8, 9, 10; Zhang 2026; de Roos Fig. 4; Bolton; Meynier; Deng 2021;
  Davidek 2006a; Xia 2022; Li 2011; Steinhagen 2021.
- **D3's pots** (Du 2023, Zhao 2019b) are M2's TEST. Their route pattern is read from the published tables; no engine
  number on them before M2's item 7.
- **Band sources** (sec. 5.2) never decide, and decisive units never set a band.
- **Not simulable:** Ames 2001 (extrusion) and Cannon 2023 (no recovery, no NF). Evidence only.

### 4.2 The decisive rows

Directions follow from each structure. "[dry]" marks arithmetic on printed values.

| # | measurement (record) | value basis; engine status | R1 | R2 | R3 (control) | R4 (control) | R5 | separates |
|---|---|---|---|---|---|---|---|---|
| D1 | **NF / MFT in the glucose pot** (Hofmann 1998 Table 5: 7.1), with the ribose pot's 2754 as in-sample context | within-pot ratio of two printed values, no tau; glucose cells **run under R0 only**, ribose cells seen (fit rows) | glucose's 1-DPO feeds NF and MFT as ribose's does, so glucose NF / MFT equals the DPO-derived ratio M2 sets. [dry] The measured NF ratio glucose / ribose is 2.5e-4 and the MFT ratio 0.096, 390x apart, so R1 needs most ribose NF to come from pentose-only routes | glucose MFT bypasses NF; NF / MFT is set only by glucose NF from other routes | NF from HA + MGO (Table 6) over MFT from HA + MP (Table 8); follows MGO's standing level, hence Kanzler's law | as R1 | **as R1** (the same 1-DPO entry and downstream shares) | R2 from R1 / R5 (scale-free; nearest the intermediate). **Does not separate R1 from R5** |
| D2 | **glucose MFT pH shape** (Hofmann Table 2: 0.3 / 1.9 / 0.4 at pH 3 / 5 / 7; ribose 55.3 / 19.8 / 2.5) | within-study contrasts log10(MFT_pH / MFT_5), no tau; pH 3 / 7 **run under R0 only**, pH 5 seen | low at pH 3 (2,3-enolisation minor in acid), falling at pH 7 as beta-cleavage takes 1-DG (Davidek 2006a) | as R1 below 5; above 5 set by the Strecker step's pH law | rising to pH 7, like fed HA + MP (Table 8) and the fragment supply (Hayashi, Yu) | set by the amine-free 1-DG supply's pH law, which no measurement fixes | set by the Amadori -> glucosone pH law, **which no measurement at >= 100 C fixes** (`k_ama_g` is licensed at pH 6.7 only); downstream as R1 | R3 from R1 / R2 by sign on the 5 -> 7 limb; the hexose shape against the pentose shape on the 3 -> 5 limb. R5 only through its declared pH sensitivity |
| D3 | **mixed isotopomers** (Du 2023, Zhao 2019b) | label share, no tau, no omega on the share (sec. 5.1); **reserved** | intact C5 | intact C5 | ~50 % mixed in a 1:1 pool | intact C5 | intact C5 | **rejects R3 by a measured share**; does not separate R1, R2, R4, R5 |
| D4 | **ribose / glucose MFT 10.4 and FFT 4.3** (Table 1), and the trend from > 870 at ~100 C (de Roos Table I) to 10.4 at 145 C | endpoint ratios; **seen, spent** | follows the entry's barrier. [dry] If the two labs' pots were comparable, the hexose route's apparent barrier would exceed the pentose route's by >= ~128 kJ/mol. They are not comparable (37 % water lactate vs 0.5 M phosphate) | as R1 | also follows MGO's loss law: with Ea 47 - 121 the C3 pool falls x2.5 - 33 faster at 145 C | the probe gave 45 - 78 at 145 C | as R1, plus R5's O2 budget | side evidence only |
| D5 | **the Cys-hexose Amadori pool** (de Roos Table II, 90 C; Du 2024b, 145 C; Zhang 2026 Fig. 3) with no hexose MFT at ~100 C (de Roos Table I, < 0.05 mg/mol Cys) | relative areas (II); a censored ratio > 870 (I); Du a C/C0 peak; de Roos I glucose arm **not run**, Zhang seen | the pool exists; its decay to 1-DG is the slow step at 100 C | as R1 | no Amadori step needed | no pool: fails structurally | the pool exists; its oxidative exit needs O2 | R1 / R2 / R5 from R4. Screen signatures (sec. 5.5): development grade after the screen |
| D6 | **FFT vs furfural vs pH** (Meynier Table 3; Hofmann F25#8) | within-study contrasts; seen | (part E, deferred) | | | | | deferred with part E |
| D7 | **amine class at fixed sugar**. Xia 2022: TriG-ARP x5 - 7, separate pots. Luo 2022: in-pot peptide / Gly-ARP ~1 at 80 - 100 C. de Roos Table II: Fru-Cys / Fru-Gly 1 -> 4 - 5 with dilution. Brands 2001: casein leaves glucose loss unchanged. Kanzler 2026: amines catalyse MGO loss without consumption. Yu 2017: within-class x2 - 7 | within-study ratios; Xia seen; Luo, de Roos II, Kanzler not run | (A-R1) each class its own constant and dilution law; catalysis declared | -- | (A-R3) equal per alpha-amine: Luo 2022 agrees at 80 - 100 C, pH 3 - 7; Xia and de Roos's dilution law disagree | -- | -- | **does not separate A-R1 from A-R3 on disk**: the same network's two studies point opposite ways, and within-class spread is as large as the contrast |

**The plan's decisive data, checked.**
- "Hofmann ribose / glucose ratios (MFT 10.4, FFT 4.3)": D4 here, seen and spent, side evidence.
- "The FFT pH slope": D6, deferred with part E.

### 4.3 The decisive data on disk, or their absence

- **R0 pinned:** D1 and D2 at refusal scale. R0's no-refit residuals there are -6 to -25 dex, so R0 cannot be favoured
  over any rival that reaches the data's scale.
- **R3 (control):** rejected by D3's measured share; D2's 5 -> 7 limb agrees in sign.
- **R4 (control):** fails D5 structurally.
- **R1 against R2 (the claim against the bypass): one decisive unit, D1.** It is unseen under every rival but R0: one
  laboratory, one temperature, NF by an unstated method, SIDA MFT.
  - Whether it decides depends on how far apart R1 and R2 put glucose NF / MFT. Under omega it needs >= 1.16 dex (x14)
    between the rivals' predictions (sec. 5.4).
  - D2 does not separate R1 from R2 (both predict the same sign on both limbs).
- **R1 against R5 (C-6 loss through 1-DG against C-1 loss through glucosone): no decisive unit on disk.**
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
- On disk, M4's part G verdict can at best be "R1 or R5 (or R2), not R0, R3 or R4". It reads "R1 or R5, not R2" only
  if the power check shows D1 separating R2 from them.
- The claim cannot be favoured over R5 without the designed experiment (sec. 4.4, arm 1 with N2 / air, or arm 2's fed
  glucosone).
- Part A stays as declared bands.

### 4.4 The designed experiment (W6)

One laboratory, three arms. Throughout:
- n = 3 independent heatings per cell;
- sealed vessels with stated fill and headspace;
- pH stated on the cooled scale for every sample, and in situ where an electrode allows;
- every assay's LOD / LOQ and recovery reported;
- SIDA with labelled standards for MFT, FFT and NF.

| arm | conditions | species | decides |
|---|---|---|---|
| **1. which carbon** | Hofmann's design: Cys 33 + glucose 100 mM, 0.5 M phosphate pH 3 / 5 / 7, 145 C, 20 min; glucose as [1-13C], [6-13C] and 1:1 [13C6] / [12C6] (and unlabelled); **each under air and under N2 (degassed, argon headspace)** (new 2026-11-06) | MFT, NF, FFT, furfural, MP, each isotopologue-resolved | see the outcome list below |
| **2. fed C6 dicarbonyls** | synthetic 1-DG (Glomb 2000 route) **or synthetic glucosone** (new 2026-11-06), 2 mM, + Na2S or Cys 6 mM, +- glycine 33 mM, 0.5 M phosphate pH 5 and 7, 145 C, 2 - 20 min, N2 and air arms | MFT, NF, DDP, 1-DPO and 3-DPO (OPD quinoxalines), acetic and formic acid, residual dicarbonyl | see the outcome list below |
| **3. fed Cys-hexose Amadori** | synthetic D-fructose-L-cysteine (de Roos's preparation) 10 mM, 0.5 M phosphate pH 3 / 5 / 7, 100 and 145 C, 0 - 60 min, N2 and air | the Amadori compound, free Cys, glucosone, 1-DG and 3-DG (OPD), MFT, NF | part G's supply constants (formation, decay, Cys return, oxidative exit) and the mechanism of D2's pH maximum |

What each outcome decides:
- **Arm 1, [1-13C] gives a labelled MFT methyl and [6-13C] unlabelled MFT:** C-6 lost, so R1 or R2.
- **Arm 1, the reverse:** C-1 lost, so R5.
- **Arm 1, ~50 % from both:** C-1 / C-6 equilibration (Tressl 1999 Fig. 14), and the N2 / air contrast decides. R5
  needs O2; R1 and R2 do not.
- **Arm 1, NF's pattern equal to MFT's:** R1 or R5. Independent of it: R2. A mixed pattern: R3.
- **Arm 2, MFT forms from fed 1-DG without the amino acid, through 1-DPO and NF:** R1.
- **Arm 2, MFT needs the amino acid, with a 1,4-dideoxy quinoxaline and no NF:** R2.
- **Arm 2, MFT and C-1-lost 1-DPO form from fed glucosone, under air or not:** R5.
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
| D1 glucose NF / MFT (Hofmann T5) | decisive (R2 vs R1 / R5) | within-pot ratio | triplicates "differ by <= 10 %" (SIDA, MFT). NF method not stated: the same 10 % as a declared transfer, sensitivity x2 | omega per printed value (2 omega^2 on the ratio); no tau | MFT, NF: none printed; declared transfer 1 % of the lowest printed value in the table, required agreement over x0.1 - x10 (sec. 5.3) | exact log-normal; censored form exact |
| D1 ribose NF / MFT | context (in-sample) | within-pot ratio | as above | as above | as above | reported, not scored |
| D2 glucose MFT pH 3 / 7 against pH 5 (Hofmann T2) | decisive (R3 by sign; the pentose-shaped null) | two within-study contrasts sharing the pH 5 point (bivariate, correlation 0.5) | <= 10 % | omega per printed value; no tau | as D1 | exact bivariate normal in log10; censored where a rival predicts below the floor |
| D3 isotopomer shares (Du 2023, Zhao 2019b) | rejects R3 | label share | Du: no SD, two blind readers identical; Zhao: triplicate means | binomial-type share error from the printed precision (+-1 % absolute) plus a mixed-channel floor of 1 %. Omega is a level error and does not apply to a share | M+1 - M+4 channels: 1 % (the printed "< 1" cells) | exact |
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
- Hofmann Tables 2 or 5 (decisive);
- Du 2023 or Zhao 2019b (reserved);
- the screen's pots (de Roos Tables I and II, Hofmann Table 1, Bolton, Du 2024b).

Starts are drawn inside the bands, never set from a decisive or screened value. Some band sources have engine outputs in
other records (Li 2011, Deng 2021, Martins 2003b, Davidek 2006a, Zhang 2026; sec. 0); none is scored or screened here. A
seen source sets a centre or a fixed measured constant, and any width it does not itself print is unsourced and carries
required agreement. First-reader sources carry the rule stated at the top.

| coordinate | band | source | kind |
|---|---|---|---|
| Glc + Cys -> Cys-hexose Amadori, formation | centred so the predicted free-cysteine first-order loss at Li 2011's conditions (0.2 M each, pH 7, 100 - 130 C) matches its printed k (0.0011 /min at 100 C), x0.3 - 3; Ea 86 - 136 kJ/mol (Zhang 2026's ~101 to Li's 121, +-15) | Li 2011 (glucose + Cys itself; free-Cys loss is a lower bound on net formation); Zhang 2026 (BTBU, units not printed: the barrier only) | same pot class; the width is **unsourced: required agreement** at both ends. Replaces the 2026-11-05 transfer from the trunk's Glc + Gly constant, which is kept as the required-agreement alternative |
| Cys-hexose Amadori -> 1-DG + Cys, decay | Deng 2021's Met-glucose ARP decay 0.0048 /min at 100 C (label pH 5.5), pH slope from its 7.5 arm, x0.3 - 3; **Ea 34 - 107 kJ/mol** | Deng 2021; Zhang 2026 (85.7); Martins 2005 (97 - 107); **Luo 2023 (Gly-ARP, 34 - 56; first reader)** | same class (hexose amino-acid ARP); the width unsourced: required agreement. Luo's decay is fast-then-slow; a two-pool form is not a rival here (no redesign) |
| Cys returned per Amadori decay | **0.54 - 1.0** | Martins 2003b, Davidek 2002 / 2005 (71 - 96 %); Deng 2021 (~1); Luo 2023 (0.54 - 1.0; first reader) | measured, same class |
| `r_odg_dpo` (1-DG -> 1-DPO + C1; R1) | 0 to a ceiling: the share of fed 1-DG not accounted for by the measured exits (Davidek 2006a; Voigt 2009; Smuda 2010) at their own conditions; barrier 40 - 135 kJ/mol | measured exits bound it from above only; Tressl 1994 and Glomb 2010 report no 1-DG C5 (no floor) | numerical above zero: **a verdict that needs a value off the floor does not ship without arm 2** |
| Cys-hexose Amadori -> glucosone + Cys (R5) | B21's shipped `k_ama_g` (Amadori -> glucosone + amine) at its value +-0.5 dex (its declared transfer LacLys -> fructosyl-glycine); pH-flat as shipped (licensed at pH 6.7), sensitivity at the 2,3-enolisation slope (Beck / Nedvidek, ~0.44 dex/pH); O2 factor 1.0 (air), sensitivity 0 once the cysteine budget closes the O2 (sec. 3) | the shipped constant; no Hamzalioglu value read here | shipped constant, transfer declared; the pH and O2 laws are **unsourced: required agreement** |
| glucosone -> formic + C5 enediol -> 3-DPO / 1-DPO (R5) | 0 to a ceiling: glucosone's total measured loss at 50 C, pH 7.4 + Lys (Glomb 2010: t1/2 8 h, k ~1.4e-3 /min); barrier 40 - 135 kJ/mol; the 3-DPO : 1-DPO split as the network's pentose enediol split; on the one glucosone pool, competing with `r_g_go` | Glomb 2010 (first reader; ratio "as in glucose reactions", no numbers) | numerical above zero: **a verdict that needs a value off the floor does not ship without arm 2's fed glucosone** |
| R2's Strecker / reduction step on 1-DG | 0 to the trunk's dicarbonyl + amino-acid Strecker class constant | M3's class. M3 phase 1 ended undecidable on disk (no verdict), and M0's MF-29 finds route A's rate law, carried from fed-dicarbonyl pots, ~5x high in the argon pot, so the upper edge is suspect | class bound; required agreement |
| NF and DDP shares downstream of 1-DPO | as M2 leaves them; if M2 has no verdict, M2's rivals carried as a factor | M2 | inherited |
| **MGO loss law (`k_mgo_mel` and any MGO sink on the network)** | **Ea 47 - 121 kJ/mol**, with the level at 100 C pinned to B60's 0.0228 /min (an amine-pot level; Pfeifer's amine-free 130 C value is x7 lower); order 1 primary, order 2 a sensitivity (identified within Kanzler 2021's series); B60's Ea 20 a required-agreement sensitivity throughout | Kanzler 2026 (MGO + amino acids, 95 - 121; second reader done for CF-122); Kanzler 2021 with Bork 2022 (MGO alone, 47 - 66; first reader); Pfeifer 2013 (130 C level check; first reader) | measured class, **one laboratory lineage**; within-lab between-paper spread 0.51 dex |
| amine catalysis of MGO loss (part A) | excess k proportional to [amine]^0.75 - 0.84 (Ala); Lys x5 at 100 mM | Kanzler 2026 | one laboratory; part A evidence; no class constant fitted on flavour rows |
| R3's HA + MGO -> NF | Hofmann Table 6 yields at pH 3 / 5 / 7 | Hofmann 1998 Table 6 (not run, not scored; the same paper as D1 / D2 but other pots), a fixed measured constant | measured constant |
| R3's fragment-supply pH law (`k_glc_ha` vs pH) | **0.28 - 0.76 dex per pH unit, rising** | Hayashi 1986 (acid limb >= 0.28, C3 0.33 - 0.41); Yu 2017 (MGO standing level 0.50 - 0.76 on 5 -> 7); both first reader | measured range, standing levels not supply constants: required agreement at both ends |
| part A class constants | peptide / amino-acid ARP formation ratio 0.43 - 7 (Luo 2022 in-pot to Xia 2022 separate pots); Lys / Gly dicarbonyl factor 2 - 7 (Yu 2017); protein lysine from B20 | Luo 2022, Xia 2022, Yu 2017, B20 | declared bands per class; never fitted on flavour rows |
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
  point makes the bivariate form a little weaker).
- **D3:** a 1 % mixed channel against R3's ~50 %: decisive at any omega (a share, not a level).

**What the synthetic power check must simulate (required before the freeze; on the first build of the network
rivals, W3).** Synthetic data only. R0 at its shipped constants may be evaluated on the designs (as the M2 prereg did).
No measured decisive value is read.
1. **Generating rivals:** R0 (pinned and refitted), R1, R2, R5, and the controls R3 and R4. Also two mixtures:
   R1 + R5 at 0.5 each, and R1 + f R3 with f = 0.05 and 0.1 at pH 7. The mixtures show whether a unit misattributes a
   mixed truth.
2. **Draws:** Monte Carlo over every sec. 5.2 band, including both ends of each required-agreement range, the first
   readers' pre-dossier bands, and:
   - M2's rivals as a factor (D1's R1 / R5 prediction depends on M2's NF / DDP shares);
   - the MGO law at Ea 47, 121 and 20, at orders 1 and 2;
   - R5's O2 factor at 1 and at the cysteine-closed budget;
   - R5's pH law flat and at the 2,3-enolisation slope.
3. **Designs, exactly as printed:**
   - D1 and D2: Hofmann's aqueous arm. Cys 3.3 mmol + glucose 10 mmol in 100 mL 0.5 M phosphate, pH 3 / 5 / 7 label,
     145 C, 20 min, 200 mL autoclave.
   - D3: Du 2023's (Cys 1 mmol + 1:1 [13C6] / [12C6] glucose, 0.2 M phosphate pH 6.5, 150 C, 120 min) and Zhao 2019b's
     GSH design, through the isotopologue pools of sec. 6.
4. **Noise:** sec. 5.1 exactly.
   - omega 0.419 per value (0.390 as a sensitivity), and the NF-assay noise at x1 and x2;
   - floors at x0.1, x1 and x10 of the declared transfer;
   - D3's +-1 % absolute and 1 % mixed floor;
   - the correlation of 0.5 between D2's two contrasts.
5. **Scoring:** the same densities, adequacy rule, pinned-R0 rule and Holm correction as item 6.
6. **Report:** for each unit and each pair, P(the generating rival wins by > 1.92 nats), and P(a "both inadequate"
   reading).
   - A unit whose power is < 0.5 for a pair is demoted to "ranks, does not decide" for that pair. This happens before
     any measured decisive value is compared with a rival.
   - The expected outcome is written in sec. 4.3: D1 has zero power for R1 vs R5 by construction.
7. **Cost:** ~4 - 6 h of work plus ~2 - 6 h of unattended compute, once the network rivals exist.

**P(the decisive data on disk can separate the rivals).**
- R0 pinned from the rest: ~0.97.
- R3 out by D3: ~0.95 (the Schieberle 2003 transfer caveat).
- R4 out by D5 (structural): ~0.9.
- R2 from R1 / R5 on D1: ~0.3. That is P(the review keeps D1 decisive) ~0.6 times P(D1's power >= 0.5 over the bands)
  ~0.45, which is low mainly because M2's shares are carried as a factor.
- R1 from R5: ~0.02.
- Part A: ~0.05.
- P(part G ends with one rival favoured over all others on disk): **~0.02**. The most likely verdict (P ~0.7) is
  development grade: "R1 or R5 (or R2), not R0, R3, R4".

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
added, promoted or redesigned around the screen. The screen's pots are then declared seen; D1 and D2 are untouched.

**Stop rule.** If no rival among R1, R2 and R5 passes S1, S2 and S5 inside its bands within ~15 h, part G reads "hexose
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
- **Amine classes as species, not a convention:**
  - free alpha-amino acid (with identity factors), cysteine, peptide (GSH, diglycine) and protein-bound lysine (one
    key, A2b);
  - each declares which sugars it competes for, so `compose` refuses undeclared competition;
  - a separate catalytic term: amine catalysis of dicarbonyl loss without consumption (Kanzler 2026).
- **One pot with all of them:** Glc + Cys + Gly + protein in one integration (today the lanes refuse to combine).
- **In-situ pH on the hexose Amadori and enolisation steps,** with the Cys-hexose Amadori compound's pKa in the one
  titratable table; the cooled-scale observable for label-pH comparisons; and a per-step hook for M5's buffer factor
  (every Hofmann pot is 0.5 M phosphate).
- **Isotopologue pools:** labelled and unlabelled copies of the sugar carbon, position-resolved at C-1, C-2 and C-6 (or
  atom tracking), so D3 and arm 1 are observables, not hand calculations.
- **The forced-input driver** (A2 contract sec. 2): 1-DG / 3-DG courses (Deng 2021) clamped while the entry is tested.

## Open questions for the build

1. **Order with M1 and M2.** Part G feeds the pentose steps whose NF / DDP shares M2 refits and whose supply M1's C1
   replaces. D1 is read on the network after M2's verdict, or with M2's rivals carried as a factor.
2. **The Hofmann glucose pH 3 / 7 arms.** They are DEV under both frozen splits, but declared Maillard-path hold-out
   bundles. This case treats them as decisive (plan sec. 7: held-out status is the split files'), scored once in item
   6, never a fit row. If the older declaration binds, D2 moves to item 7 and R2 vs R1 / R5 rests on D1 alone.
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
9. **D1 and D2 after R0's probe.** Both pots were run under R0 (the no-refit thiol sub-model probe) before this case,
   and the authors read that probe's per-data-set summary while refreshing. This case keeps them decisive between R1,
   R2, R3, R4 and R5, which have never run there. If the review rules that any engine number on a pot spends it, part
   G has no unseen decisive unit on disk. Its verdict is then development grade until the designed experiment or the
   sweep-2 fed-ARP candidate supplies one.

## What should block the freeze (for the review and the coordinator)

1. **The synthetic power check cannot run yet.** It needs the network rivals, and the network infrastructure (PR #74)
   and the sulfur port are not merged. The freeze waits on W3, by the plan's own dependency (W3 before M4).
2. **Open question 9 (does R0's probe spend D1 / D2?).** This is a ruling, not a computation. If yes, part G has no
   decisive unit on disk at all.
3. **Open question 2 (the hold-out declaration of the glucose pH 3 / 7 bundles).** This is a ruling. If the older
   declaration binds, D2 leaves the decisive set.
4. **R5 is new.** The review should test whether R5 is a fair rival or a straw man. Points to test: whether glucosone
   from the Cys-hexose Amadori compound at 145 C is plausible with cysteine scavenging the O2, and whether B21's
   `k_ama_g` (licensed at pH 6.7, fitted on LacLys) may be used at its shipped value. It should also test whether the
   claim should be restated as "a hexose C5 entry into 1-DPO, carbon origin undecided" (R1 + R5 as one claim, with the
   designed experiment to split them). That restatement would change item 1, so it must be settled before the freeze.
5. **Not blockers, but the owner should know:** on disk, the most likely verdict is development grade (P ~0.7). The
   claim cannot be favoured over R5 without W6, which comes last by the owner's ruling (2026-11-05).

## Proposed update to `docs/reference/decisive_measurements.yml` (entry `glucose_mft_carbon_origin`)

To be applied at the freeze by that file's one writer; not applied here.
- Replace the rival "unlisted: a route through a C-1-lost C5 (no source on disk proposes it ...)" with "R5: the
  oxidative route, Cys-hexose Amadori -> glucosone -> beta-cleavage (C-1 lost as formic acid) -> 1-DPO -> the pentose
  steps (Glomb 2010)".
- Mark R3 and R4 as declared-failing controls.
- Arm 1: add the air / N2 sub-arms (o2: "air and N2 arms, stated headspace").
- Arm 2: add synthetic glucosone beside 1-DG, with 3-DPO and 1-DPO measured.
- Outcomes: "C-1 lost" now means R5, not "an unlisted 3-DG-side route". "~50 % from both labels" means equilibration,
  and the N2 / air contrast decides.
- `why_not_on_disk`: add "Glomb 2010 (Halle, 50 C) labels the hexose 1-DPO as C-1-lost via glucosone; no source at
  >= 100 C measures it".

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
