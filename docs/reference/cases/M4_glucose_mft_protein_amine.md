# Case M4: glucose to MFT; protein and other amines as classes

**REFRESHED 2026-11-05; items 1 - 5 freeze after the sweep's candidates are read or excluded** (sec. 2.9) **and after
the power check of sec. 5.4.** An independent review follows: one round plus one confirmation. Nothing here is frozen.

*Plan step M4 ([STATUS_AND_PLAN sec. 5, W2](../../guides/STATUS_AND_PLAN.md); errors E1 (F23-a, F24-a) and E8 (F25-a,
-b, F06-b) of [ENGINE_VS_LITERATURE](../pathway_atlas/ENGINE_VS_LITERATURE.md) sec. 3, amended in sec. 5 (E1: "engine
error SUPPORTED; band-now fix REFUTED"), and the plan item "amine as a variable" (F02-b, F05-a, F19-e). Items 1 - 5 of
sec. 7, first written 2026-11-04 and refreshed 2026-11-05 to this week's rules (M1.1's measurement model, the item-2
citation sweep, the M1c template), before any result of the change exists. No engine was run for this case, on any pot.
Items 6 - 9 are appended later and never edited into this part. M4 is built on the one network (W3,
[A2_MODULE_CONTRACT](../../guides/A2_MODULE_CONTRACT.md)); sec. 6 says what it needs from it. Part E (E8) is kept for the
record but is not built: the plan's M4 row defers it until M1's cysteine-Amadori route is decided.*

**Rules this case applies (M1.1, 2026-11-05).**
- **Value basis per series.** Each decisive series states whether it is an absolute level, a within-pot ratio, a
  within-study contrast, a C/C0 course or a label share (sec. 5.1). H0's between-laboratory tau (tau_ABS 0.964 dex) enters
  only on absolute levels compared across laboratories or with an engine level, never on a ratio, contrast, C/C0 or
  share.
- **omega = 0.419 dex** (= 0.965 ln) per printed value in every unit's primary covariance; 0.390 is a sensitivity.
- **Floors** come from each assay's reported LOD / LOQ or reporting limit, never from a data point. Where a paper prints
  none, the floor is a declared transfer with required agreement (sec. 5.3).
- **Densities** are exact, or Monte-Carlo with PSIS k-hat < 0.7 and ESS >= 100, on every form that can flip a sign,
  including censored (below-floor) forms and every sensitivity, not only the primary.
- **Adequacy.** A unit favours a rival only where that rival is adequate (posterior-predictive p > 0.05 under the
  unit's own measurement model) or the rivals' means differ by more than omega; otherwise the unit reads "both
  inadequate". Holm across guard groups.
- **The pinned incumbent.** A favoured rival must also beat R0 pinned at its shipped (B60) constants, not only R0
  refitted.
- **Bands and starts** never come from data scored or screened in this case (sec. 4.1, 5.2). An unsourced number enters
  with required agreement: the verdict must hold at both ends of its declared range.
- **Screen first, demotions only.** Before any fit, a mechanism-signature screen (sec. 5.5); after it a rival can be
  demoted, never promoted or redesigned, and the screen's pots are declared seen.
- **Power check** before the freeze (sec. 5.4).
- **pH scale.** Unless stated, every pH below is a **label pH** (the buffer or solution as prepared, at room
  temperature, before heating). "Held" = pH-stat or a buffer with a reported drift; "cooled" = read after heating and
  cooling; "in situ" = at temperature (an engine quantity unless a paper says otherwise). The engine is compared with a
  measurement only on the measurement's scale (`observations.ph_cooled` for label or cooled readings).

**Split notation.** "v1" = [split.yml](../contrast_benchmark/split.yml); "v2" =
[split_v2_per_change.yml](../contrast_benchmark/split_v2_per_change.yml). v2 E1: DEV labs BTBU, CAAS, DFA Garching,
Nestle / Cerny; **no TEST contrast** (UNTESTABLE_on_held_out). v2 E8: DEV labs DFA, Reading, Wageningen; **no TEST
contrast**; Cerny 2007's F13#7, F25#6, F25#7, F25#22 excluded as calibration input. Atlas cells are 0-based indices of
the F files in [pathway_atlas](../pathway_atlas/). Every cell below is v1 DEVELOPMENT unless marked. Two cautions:
- F23#7, F23#8, F24#4, F24#5 and F25#23 (Zhao 2019b, Du 2023) are DEV(E1) but **TEST for E6 and E11 in v2** (M2's
  errors), and their keyed tables (`heldout_du2023_tab2.yml`, `heldout_zhao2019b_tab1_2.yml`) are TEST-only
  ([contrast-benchmark README](../contrast_benchmark/README.md) sec. 11). This case reads only their route pattern
  (intact C5 vs mixed isotopomers); anything they say about the norfuraneol share or the route mix with temperature is
  M2's held-out evidence and is read in item 7 only. **No engine run on them is allowed before M2's item 7.**
- The Hofmann glucose pH 3 / 7 arms are DEV under both frozen splits, and separately declared Maillard-path hold-out
  bundles (`data/benchmarks/external_validation/maillard_path/mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH{3,7}.json`);
  the pH 5 arm is a trust-loop predictive row, not a B60 fit row (`kinetic_core_b60_fit_targets.json` holds only the
  ribose and xylose T1 rows). See open question 2.

## 0. What has been seen under the engine

**By the author.** Engine outputs read:
- ENGINE_VS_LITERATURE Appendix A probes C and E (glucose + Cys MFT 1.5e-14 mM at pH 6, 100 % C2 + C3; ribose / glucose
  MFT 4.4e10 and FFT 0.47 on the Hofmann T1 design, pH 5) and B / D (ribose + Cys FFT and MFT vs pH 3 - 7, 30 mM);
- the sec. 5 red-team probe of `HEXOSE_C5_ENTRY` (ribose / glucose MFT 45 - 78 without amine, 0.79 with 33 mM glycine,
  T1 design; a scratchpad probe, not an artifact);
- `core_panel_scores.md` rows for the three Hofmann glucose bundles (all "refused", no number) and
  `thiamine_cys_glucose_120C_Bolton1994` (predicted 232.9 vs 11.7 ppb);
- **read while refreshing (2026-11-05):** the per-data-set summary of the no-refit thiol sub-model probe
  (`kinetic_core_thiol_submodel_probe.md`, shipped constants = R0): Hofmann Table 5 NF mean residual -12.3 dex over 2
  rows (worst -24.8), Table 2 MFT -6.0 over 4 rows (worst -9.7), Table 1 MFT -5.4 over 4. So D1's and D2's pots have
  been run under R0 and that summary is seen; no per-row number was read, and no rival other than R0 has run on them.

No engine output exists on any v2 E1 or E8 TEST contrast (both sets are empty). Measured values only for the v1 TEST
cells F06#16 (Lee 2022), F14#5 (Hollnagel 2000) and the v2 E5 TEST cell F06#11 (Thornalley 1999), read in the atlas
files. The sweep (sec. 2.9) was screened on titles and abstracts only.

**Pots by engine status.** "Seen" = an engine number on that pot exists in a record (fit row, probe or screen), so the
pot cannot be decisive here. "Run under R0 only" = the incumbent's no-refit probe has a number there and no other
structure has run; such a pot stays decisive between the new rivals (R0's numbers are refusal-scale and predate every
rival's design), which the review is asked to check (open question 9).

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
| Bolton 1993 | **seen** | core panel |
| Deng 2021; Davidek 2006a; Martins 2003b, 2005; Brands 2001, 2002c | **seen** | pH-acid sub-model probes; trunk fit rows |
| Xia 2022 | **seen** | B21 |
| Tressl 1994; Voigt 2009, 2010; Smuda 2010; Kanzler 2026; Beck 1988; Nedvidek 1992 | **not run** (no engine species for most of their observables) | -- |

## 1. The claim

Three parts, built and judged separately on the one network.

| part | step(s) | today (B60, read from the code) | M4 |
|---|---|---|---|
| G. hexose entry to MFT | new: Glc + Cys -> cysteine-hexose Amadori compound (D-fructose-L-cysteine) -> 1-DG + Cys released; then `r_odg_dpo` (1-DG -> 1-DPO + C1, [sulfur.py](../../../src/kinetic_core/sulfur.py) ~927, "the C-6 loss") into the shared pentose MFT steps | glucose reaches MFT only through `r_glc_c2c3` (log10 k -8.59 at 145 C; Glc -> HA + MGO) then `r_mgo_mp` (MGO + H2S -> 1-mercapto-2-propanone, MP) and `r_ha_mp_mft`; FFT through `r_glc_fur` (`acid` factor); `HEXOSE_C5_ENTRY.k_per_min_100c = 0` ([parameters_sulfur.py](../../../src/kinetic_core/parameters_sulfur.py) ~623 - 641); cysteine is not an amine on the trunk, so a Glc + Cys pot has no amine route to 1-DG; `engine.unidentified_routes` refuses hexose-only thiol requests (engine.py ~1024 - 1043) | a measured-class Cys-hexose Amadori formation and decay, `r_odg_dpo` on with a declared barrier |
| E. the furfural supply's pH law (E8) | `r_pent_tdp`, `r_arp_tdp`, `r_glc_fur` (`acid` = [H+], a decade per unit; sulfur.py ~1160 - 1175) | as stated | an amine-catalysed 1,2-enolisation term, pH-flat or weakly rising, beside a weaker acid limb. **Deferred** (plan M4 row) |
| A. amine classes | the trunk's `r_schiff` (Glc + Gly), the sulfur lane's Cys-catalysed `r_pent_dpo` / `r_pent_tdp`, B20's `r_glc_lysp` (protein lysine, trunk only) | one amine: MET and PRO folded into `Gly` (engine.py ~2560 - 2570); `LYSP` exists only on the trunk lane; no peptide | classes: free alpha-amino acid, cysteine, peptide, protein-bound lysine; each with its own condensation constant; their competition for each sugar declared, and **amine catalysis of dicarbonyl loss declared separately from amine consumption** (Kanzler 2026, sec. 2.2) |

**Part E depends on M1.** Enabling `CYS_AMADORI_ROUTE` (M1's combination C1) retires `r_pent_dpo` / `r_pent_tdp`
(`inactive_reactions`, sulfur.py ~975 - 978) and carries the pentose supply on de Roos 2005 Fig. 4's measured,
U-shaped Cys-Amadori decay. If C1 ships, E8's pentose part is decided by M1 / M1b, and part E keeps only
`r_arp_tdp` (fed ARP pots) and `r_glc_fur`.

**MGO on the one network (new 2026-11-05).** Today the sulfur lane has no amine-free MGO sink: its MGO (from
`r_glc_c2c3`, `r_dpo_c2c3`, `r_nf_decay_oh_mgo`, `r_tdp_mgo`) leaves only by `r_mgo_mp`, `r_cys_actz` and the
switched-off `r_mgo_cys_*`. On the one network MGO is one species, so the trunk's `r_mgo_mel` (MGO -> melanoidin,
first order, amine-free, **fitted**: log10 k(100 C) -1.643 +- 0.114, Ea **20.0 kJ/mol at its search floor**,
`FITTED_EA_BOUNDS` (20, 260)) acts on every MGO the sulfur module makes. Kanzler 2026 (TU Berlin, sec. 2.2) measures
MGO loss at pH 5 (initial, unbuffered) with an apparent Ea of **95 - 121 kJ/mol** (MGO + Ala / Phe / Pro, 80 -> 100 C;
~95 - 113 kJ/mol for MGO alone against Bork 2022, same laboratory). From the same 100 C value, an Ea of 95 - 121 gives a
loss x13 - 33 faster at 145 C than Ea 20 does [dry]. Every rival that carries MGO (R3 most of all; R1 and R2 through
NF's MGO-releasing decay) therefore inherits a sink whose temperature law the only measurement contradicts. Sec. 5.2
bands it; R0 is also scored pinned at Ea 20.

**Earlier waves on these structures** (cited from the record, none re-run):

| structure | wave | outcome | kind |
|---|---|---|---|
| hexose MFT through `k_glc_ha` (C2 + C3), carried by the eight Hofmann T1 level rows | B8 | the route existed only through the rows scored on it | whole-model fit |
| the same with the T1 level rows dropped | B9 (`tasks/data_restructure_plan.md`, 2026-09-03) | glucose and fructose MFT predict zero; `k_glc_ha` at its band floor; B9's Laplace covariance singular in that coordinate | whole-model fit |
| the hexose route reopened | B47, B51.1 | DO NOT SHIP: `k_glc_ha` +3.6 decades bought back one Whitfield row; directional claim SUG-09 flips | whole-model fit |
| hexose -> C5 entry (`HEXOSE_C5_ENTRY`) | sub-model 4 (backlog, unchecked) | never pre-registered or run | -- |
| the entry switched on alone (the "band-now" fix) | red team, ENGINE_VS_LITERATURE sec. 5 | ribose / glucose MFT 45 - 78 with no amine (measured 10.4), 0.79 with 33 mM glycine | read-only probe |
| protein-bound lysine glycation | B20 (SHIPS, trunk only) | Nguyen 2016's constants; answers only with `protein_g_per_l` | fit on printed rates |
| pentose + glycine entry | EX-012 / EX-013 | admitted on audit; norfuraneol sanity bound FAILED (25 - 36 mol %) | fit |
| the furfural supply's pH law | Appendix A probe B | FFT x44 from pH 4.5 to 6 and ~x400 to 6.5 against x6.1 measured (Meynier, held) | probe |
| the trunk's MGO sink | B1 (`k_mgo_mel`) | Ea at its 20 kJ/mol floor; Kanzler 2026 measures 95 - 121 | fit; conflict not yet registered |

**Errors stack (CF-107).** The B8 / B9 / B47 / B51.1 outcomes are whole-model fits and weak evidence about the hexose
route. Candidate compensating errors:
- E6 (MFT through free norfuraneol, 99.9 %) and E11 (one 64.08 kJ/mol formation barrier);
- the frozen pot-independent `k_mft_decay` and the Strecker sulfur sent to an inert pool (SULFUR_FLUX_MAP sec. 6);
- MFT destroyed in excess H2S at B9's optimum (backlog, "next sulfur question");
- E2 (the enolisation split frozen in pH) and the trunk's MGO level (B21), now with `k_mgo_mel`'s Ea at its floor
  against Kanzler's 95 - 121 kJ/mol.

The probe's 45 - 78 is an endpoint ratio that inherits E6 and E11, so it does not reject the entry either. Only the
labelling (sec. 2.1) is a measured share that makes a route inert: it rejects C2 + C3 as the hexose route.

## 2. All the evidence, by kind

### 2.1 Isotope labelling

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Du 2023 (F23#8, F24#5; ledger `du2023_table2_camola_isotopomers`) | Cys 1 mmol + 1:1 [13C6]/[12C6]-glucose, 0.2 M phosphate pH 6.5, 150 C, 120 min: MFT M+0 : M+5 = 52 : 48 (57 : 43); no mixed channel printed | hexose MFT is an intact glucose C5; C2 + C3 absent on the label basis | DEV(E1); TEST(E6, E11) |
| Zhao 2019b (F23#7, F24#4, F25#23; `zhao2019b_table2_camola_glucose_labelling`) | GSH + glucose 60 / 60 mM, phosphate pH 6.5 (molarity not stated), 140 C, 5 h: MFT ~100 % glucose C5; FFT ~94 % intact C5 | the same in a peptide-amine pot | DEV(E1); TEST(E6, E11) |
| **Tressl 1994** (TU Berlin, Tressl group; `tressl1994_extraction.md`; new here) | **positional** [1-13C]- and [6-13C]-glucose + Cys, equimolar, pH 5, 160 C, 1.5 h: every glucose-derived C5 sulfur heterocycle (furfuryl alcohol, 2-hydroxymethylthiophene, the thiolanones) is glucose C2 - C6 (C-1 lost, via 3-DG); every 1-DG product keeps all six carbons. **MFT and NF not reported** | which carbon a hexose C5 loses. The only C5 cleavage it documents (C-1 from 3-DG) cannot give MFT's 2-methyl group; no 1-DG-derived C5 (`r_odg_dpo`'s C-6 loss) is reported (absence of a report, no LOD) | DEV; a fourth laboratory |
| Liu 2023b (F23#6, F24#6; CF-086) | ribose 168 C: fragments ~6 % after the [M-H]+ correction | fragments stay minor at high T (pentose) | DEV; TEST(E6, E11) |
| Cerny 2003 / 2004, 2007b (F23#0 - #3, F24#0 - #3) | pentose MFT intact C5; free NF 7 % (pot D) | the downstream steps M4 reuses | DEV |
| Bolton 1993 (F27#6; `bolton1993_tab2_mft_34s_split`) | 34S-Cys + thiamine in a glucose bouillon, a_w 0.83, 120 C: Cys supplies ~8 % of MFT sulfur; no MFT without thiamine (< ~2.3 % of the thiamine pots) | an upper bound on Cys + glucose MFT at 120 C in a matrix | DEV |
| Davidek 2005, 2006a (F07#5 - #7) | fed DFG and 1-DG: acetic acid ~70 % from glucose C1 / C2 | 1-DG's own exit (beta-cleavage), the competitor of any 1-DG -> C5 step | DEV |
| Voigt 2010 (Halle, Glomb group; new here) | [13C]glucose (C-1, 2, 3, 5, 6) + Ala, 50 C, 7 d, air: lactic acid 48 % C1 - C3, 36 % C3 - C5, 8 + 8 % C4 - C6 | 1-DG's fragmentation pattern at low T | DEV |
| none | **no positional label on MFT itself** (which carbon glucose MFT loses: C-1 or C-6), and no time-resolved hexose label | R1 / R2 (C-6 lost, methyl from C-1) against any route through a C-1-lost C5 | gap; the designed experiment (sec. 4.4, arm 1) is Tressl's design with MFT and NF measured |

### 2.2 Fed intermediates

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Hofmann 1998 Table 8 (F24#7, F24#8) | fed HA + **1-mercapto-2-propanone** (MP; not C2 + C3 + H2S, corrected 2026-11-05), 145 C: MFT 0.01 / 0.23 / 0.27 mol % at pH 3 / 5 / 7 (rising with pH) | the C2 + C3 route's last step and its pH shape (rival R3) | DEV; B60 fit rows (`hofmann1998_c2c3_recombination_*`) |
| Hofmann 1998 Table 6 (new here) | fed HA + MGO (1 mmol each, no sulfur), 145 C: **NF 0.13 / 0.78 / 3.18 mol %** and furfural 0.11 / 0.40 / 1.60 mol % at pH 3 / 5 / 7 | **C2 + C3 make NF too**, rising with pH: R3 is not "no NF" (sec. 3) | DEV; not run |
| Hofmann 1998 Table 7 | MGO + H2S, pH 5: MP 1.8 / 4.0 mol % at 1 / 2 mmol H2S (super-linear) | the C3 sulfur step | DEV; not run |
| Hofmann 1998 Tables 3 / 4 (F13#2, F23#12, F25#4) | fed NF / furfural + H2S, pH 5 | downstream steps shared with ribose | DEV; fit rows |
| Deng 2021 (F03#11 - #12, F06#6, F07#11) | fed Met-glucose ARP, 100 C, unbuffered (pH measured only in the buffer arm): k 0.0048 (label pH 5.5), 0.012 /min (label pH 7.5); Met returned ~mole for mole | a hexose amino-acid ARP decay and its pH slope (the class analogue for the Cys-hexose ARP) | DEV |
| Martins 2003b (F03#3 - #5, F06#5, F07#10); Davidek 2002, 2005 (F03#0 - #2, F02#11) | fed DFG: amine returned 71 - 96 %; split moves with pH | the hexose ARP steps the claim reuses | DEV; Martins a calibration input |
| Davidek 2006a (F07#5, #9) | fed amine-free 1-DG, 0.2 M phosphate: acetic 0.20 / 0.56 / 0.80 mol/mol at pH 5 / 6 / >= 6.5; 1-DG gone within 5 min at pH >= 6 | 1-DG's lifetime and its main exit | DEV |
| Voigt 2009, Smuda 2010 (Halle; new here) | fed 1-DG + Ala or Nalpha-BOC-Lys, 42 mM each, 0.1 M phosphate pH 7.4, 37 / 50 C, air vs deaerated: exits are acetic, glyceric, lactic, formic acids (310 : 88 : 58 : 27 mmol/mol, 7 d, BOC-Lys), C4 dicarbonyls and their amides; negligible MGO; **no C5 product reported**; no amine-free arm | 1-DG's measured exits at low T; a ceiling on `r_odg_dpo` only by the unaccounted share (sec. 5.2) | DEV |
| **Kanzler 2026** (TU Berlin, Kanzler group; `kanzler2026_extraction.md`, first reader only; new here) | fed MGO 100 mM alone or 1:1 with ten amino acids (no Gly, no Cys), water, initial pH 5.0 (unbuffered, room T), 100 C, 0 - 120 min; 60 - 140 C at 30 min: MGO alone 65 / 31 % left at 30 / 120 min (second order fits); amino acids speed MGO loss x2.6 - 9.7 at 30 min **while barely consumed** (Lys 98 % left, Ala ~100 %); Ala order 0.75 - 0.84; apparent Ea 95 - 121 kJ/mol | (i) the temperature law of the MGO sink every rival inherits (`k_mgo_mel` Ea 20 at its floor); (ii) part A: amines catalyse dicarbonyl loss without being consumed, so a class model that only consumes amine stoichiometrically is wrong for MGO; (iii) R3's MGO standing level at 145 C | DEV (no split cell); not representable as a pot (no Gly / Cys); used as a band source (sec. 5.2) after its second reader |

### 2.3 Trapping

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Beck 1988, Nedvidek 1992 (F06#0 - #2, F07#0 - #3) | in-situ OPD on hexose ARPs / glucose + beta-Ala: 1-DG : 3-DG rises x12.5 from pH 4.5 to 7; a secondary-amine ARP goes to 1-DG even at pH 4.5 | the hexose 2,3-enolisation's pH law; amine structure moves the split | DEV (one school) |
| Nedvidek 1992 (F23#16) | the 1,4-dideoxypentosone quinoxaline needs alpha-amino acids | a 1,4-dideoxy C5 made by amino-acid Strecker on a deoxyosone (rival R2) | DEV |
| Gobert 2009 (F06#3, F07#4) | 50 C: air lowers 1-DG x3.4 | O2 on the 1-DG supply (low T only) | DEV |
| Hollnagel 2000 (F14#5) | dry glucose + OPD: deoxyhexosuloses < 1 vs 11 mol % with glycine | amine catalysis of enolisation | **v1 TEST**: listed, not used |
| Thornalley 1999 (F06#11) | BOC-Lys raises 3-DG formation x6 at 37 C | amine catalysis | **v2 TEST(E5)**: listed, not used |

### 2.4 Measured intermediates and lifetimes

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Hofmann 1998 Table 5 (`hofmann1998_t5_fa_nf_pairs`; dossier sec. 2.5) | Cys 33 + sugar 100 mM, 0.5 M phosphate pH 5, 145 C: NF / MFT 2754 (ribose), 7.1 (glucose), 23.9 (rhamnose); furfural / FFT 5.58 / 2.46 / 7.75. NF and furfural method not stated (MFT, FFT by SIDA) | whether glucose MFT passes through the intermediate ribose MFT does (D1) | DEV; ribose NF / furfural are fit rows; the glucose cells run under R0 only (sec. 0) |
| Hofmann 1998 Table 9 (`hofmann1998_t9_m2p_glucose_cys`) | 1-mercapto-2-propanone from glucose + Cys 11.4 / 59.6 / 26.5 ug at pH 3 / 5 / 7 | a C3 marker in the same pots, maximum at pH 5 | DEV; **seen** (`thiol_submodel.yml`) |
| de Roos 2005 Table II (`deroos2005_table2_cysamadori_glucose_aa`) | glucose + Cys + Gly or Pro, 90 C: D-fructose-L-cysteine forms readily; Fru-Cys / Fru-Gly 0.9 - 1.1 (15 % water), 4 - 5 (79 % water) | the Cys-hexose Amadori compound exists in situ; amine competition for glucose | DEV (Givaudan); relative RI areas |
| de Roos 2005 Table I (F23#11; `deroos2005_table1_thiol_yields_cys_carbonyl`) | 1 : 1 Cys + sugar, lactate pH 5, 37 % water, boiling 1.5 h: MFT ribose 43.5, xylose 36, glucose < 0.05 mg/mol Cys (the paper's reporting limit) | at ~100 C the hexose route is shut although its Amadori compound forms (Table II) | DEV |
| Zhang 2026 (`zhang2026_fig3_gsh_glucose_nonvolatiles`, `zhang2026_tables6_s8_volatiles`) | GSH + glucose 60 / 60 mM, 0.2 M phosphate pH 6.5, 130 - 150 C: Cys-, GSH-, Gly-Amadori by LC-MS/MS (figure-only); MFT only at 90 / 120 min (0.33 - 2.32 ng/g, RF 1) | the Cys-hexose Amadori pool in a peptide pot; MFT's late onset | DEV (BTBU); **seen** (thiol sub-model rows) |
| Hofmann 1999 (F06#4, F07#20) | xylose osones peak at 10 min, glucose at 40 min | the hexose supply's time scale | DEV |
| Brands 2001 (F14#0 - #1, F02#13, F05#0) | 3 % caseinate + glucose, 120 C: net glucose loss identical with and without casein; formic x1.4, acetic x2.6 with casein | protein-bound lysine changes products, not sugar loss | DEV |
| Xia 2022 (F02#4 - #5) | Glc + Gly / diglycine / triglycine: TriG-ARP forms 6.9x (100 C), 5.0x (130 C) faster | peptide amine as its own class | DEV (Ho network) |
| Du 2025 (BTBU; new here) | GSH + glucose 200 / 200 mM, 0.2 M phosphate pH 6.5, 145 C: MFT 135 ug/L at 120 min (5.9e-6 mol/mol GSH); CAMOLA on five sulfur volatiles, **not on MFT** | MFT level in a peptide pot (part A, weak); its lipid arms are M5's | DEV |

### 2.5 Mass balance

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Brands 2001 (F05#1, F02#12) | fructose + casein: lysine lost without Heyns or furosine | the protein amine's fate is not the free-amine pattern | DEV |
| Cao 2017 (F04#17, F04#18) | the Cys-bearing pool conserved over 5 h; glycine 16 - 24x Cys drains TTCA to ~0.03 M | amine competition in a cysteine pot | DEV (BTBU); seen (M1c) |
| Zhang 2025b (F04#19) | lysine 16x Cys moves TTCA and Cys-Amadori | the same with a free basic amine | DEV; seen (M1c) |
| Kanzler 2026 Fig. 1F | amino acids left after 120 min with MGO: Ala 100, Lys 98, Glu 90, Phe 63, Arg 47.5, Gln 15 % | which amine classes are consumed by a dicarbonyl and which only catalyse its loss | DEV |

### 2.6 One factor varied inside one study

| source (record) | factor -> result | bears on | split |
|---|---|---|---|
| Hofmann 1998 Table 1 (F23#21, F25#10; `hofmann1998_t1_cys_carbohydrates`) | ribose / xylose / fructose / glucose / rhamnose + Cys, pH 5: MFT 19.8 / 14.3 / 2.5 / 1.9 / 0.8 ug; FFT 12.1 / 9.6 / 3.2 / 2.8 / 0.8 | ribose / glucose MFT 10.4, FFT 4.3 (D4) | DEV; **seen** |
| Hofmann 1998 Table 2 (F23#18, F25#8; dossier sec. 2.2) | 0.5 M phosphate pH 3 / 5 / 7: glucose MFT 0.3 / 1.9 / 0.4, FFT 0.7 / 2.8 / 0.6; ribose MFT 55.3 / 19.8 / 2.5, FFT 22.9 / 12.1 / 1.2 | the hexose route's pH shape against the pentose one (D2) | DEV; glucose pH 3 / 7 are Maillard-path hold-out bundles, **run under R0 only** |
| Meynier 1995 Table 3 (F25#9, F13#5; CF-001 conditions) | Cys 50 + ribose 35 mM, 0.2 M pyrophosphate, held (cooled drift <= 0.2), 140 C, 1 h, pH 4.5 -> 6.5: FFT 8.5 -> 1.4 (x6.1), furfural 3.9 -> tr (> 39x), MFT 15.2 -> tr | E8: FFT falls far less than furfural (D6) | DEV(E8) |
| Hofmann 1998 Table 2 (F25#8) | ribose FFT x19 from pH 3 to 7, highest at pH 3 | E8 | DEV(E8) |
| Cerny 2007 (F25#7, F13#7) | FFT x2.3 down from pH 4 to 6; zero at pH 7 | E8, third lab | v1 DEV; v2 excluded (calibration input) |
| Schirle-Keller 1992 (Minnesota; new here) | glucose 0.5 M + Cys 0.1 M, unbuffered, 80 - 150 C: net furfural rises ~30x from 100 to 120 C; printed apparent Ea 147 kJ/mol (refit 156); no pH, no sulfur volatiles | `r_glc_fur`'s temperature law (part E; deferred) | DEV |
| Mottram 2002 (F23#20, F14#8) | buffered vs unbuffered: MFT x2.2 - 2.6, FFT unchanged; ribose-only furfural unchanged | the branches respond differently (shared with M5) | DEV |
| Blank 1997 (F13#11, F14#10) | HDMF from xylose needs glycine or alanine (< 0.01 vs 2.6 ug/mmol) | amine catalysis at 90 C | DEV; calibration input |
| Tang 2024 (F05#16 - #17) | soy vs pea isolate, wet 85 C: lysine -22 vs -8 %; the order flips dry | protein identity as a class constant | DEV (Wageningen) |
| Zhou 2023c, Feng 2022 (F27#24 - #26) | GSH a weaker MFT source than Cys; the ribose-GSH ARP makes MFT itself | peptide vs Cys class | DEV; Feng a calibration input |
| Zhang & Ho 1991b (Rutgers; new here) | glucose + Cys vs glucose + GSH, 50 / 50 mM, water, label pH 7.5, 180 C, 1 h: total volatiles 969 vs 404 mg/mol; **no MFT or FFT detected in either pot** (no LOD) | part A (Cys vs peptide at fixed sugar); a non-detect of glucose MFT at pH 7.5, 180 C | DEV |
| Martins 2010 (Unilever; new here) | 32-run D-optimal design, 100 C, 0.2 M pyrophosphate pH 6: xylose, glucose, 5-ketogluconate, Cys, Gly, Glu; MFT, FFT, NF responses as regression signs only | part A competition (Gly vs Cys at fixed sugars), weakly; no per-run data | DEV |
| Bolton 1993 (CF-090) | + 6.3 mM xylose leaves MFT unchanged in the glucose matrix | thiamine dominance in hexose foods | DEV |
| Lee 2022 (F06#16), Kocadagli 2016 (F06#15, F07#18) | amino acid raises osones 4 - 10x; NaCl moves the split | amine catalysis | **v1 TEST**: listed, not used |
| Xu 2007 (`heldout_xu2007_tab1.yml`) | pH 5.6 / 7.0 / 8.0 control arms | E8 | keyed held-out (not in a frozen split): item 7 only |

### 2.7 Fitted time courses and fits

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Zhang 2026 Table 2 (`zhang2026_table2_rate_constants`) | 19 k at 130 / 140 / 150 C; Cys + glucose -> Cys-Amadori Ea ~101, Cys-Amadori -> alpha-dicarbonyls 85.7 kJ/mol (both already B59b's transferred barriers); units not printed, no CI | barriers for part G | DEV (BTBU) |
| Martins 2005, 2005b (F06#12 - #13, F07#15 - #16) | k4 +0.35, k7 +0.60 dec / pH; Ea 97.1 / 107.3 | the hexose split | DEV; calibration input |
| Brands 2002c, Nguyen 2016, Wolf 1977a (F05#3 - #8; CF-072, CF-096, CF-100) | glycation 0.39 dex / pH; Ea 97 - 138 kJ/mol; lysine return disputed | the protein-amine class constants | DEV |
| engine records | B8, B9, B47, B51.1, B20, B1 (`k_mgo_mel`), EX-012, EX-013, red-team probe | sec. 1 table | engine records |

### 2.8 Left out, with reasons

- The thiamine route (F27, E18), except Bolton's bound: a separate MFT source with its own case.
- Acrylamide amine competition (Claeys 2005, Rydberg 2003; F29#3, #21): E9, the acrylamide lane.
- Pyrazine amine identity (Leahy Lys >> Asn, F19#10; Van Lancker dipeptides): M6.
- Strecker aldehydes from glucose pots: M3. Lipid with glucose (Du 2025's fatty-acid arms, Zhao 2019b Table 1, Yang
  2015): M5.
- Goncuoglu 2026 PPI extrusion (F29#24, v2 TEST E9) and Hamzalioglu 2026 (F05#10 - #11, v1 TEST): protein amine in
  acrylamide / milk chemistry; held out.
- Hofmann 1997 (glucose / rhamnose + Cys, the Table 1 pots): FD factors only, no amounts. Hofmann 1996b (Weurman):
  declared duplicate of Hofmann 1998's tables (`pdf_dossier_aliases.yml`).
- Madruga 1995 (beef psoas, pH 4 - 5.6, MFT rises in acid): a meat matrix whose sugar is mostly IMP-derived ribose, not
  a hexose route; judgement: W5 (recipe) evidence, not M4's.
- Liu 2021 (FFT from glucose + Cys with ribose): part E, deferred.

**Counts by kind (sources).** Labelling 8 (one gap); fed intermediates 9; trapping 5 (two held out); measured
intermediates 10; mass balance 4; one factor within a study 16 (three held out); fitted 4 plus the engine records.

### 2.9 The citation sweep (sec. 7 item 2), 2026-11-05

**Keyed to M4's decisive criterion, not its topic.** The question asked of each work: does it hold (a) NF and MFT
measured in one hexose + cysteine (or H2S) pot (D1); (b) hexose + Cys MFT at >= 2 pH in one study (D2); (c) a
positional glucose label on MFT, or a fed hexose intermediate (1-DG, the Cys-hexose Amadori compound) giving MFT; (d) an
amine-class contrast at fixed sugar (D7); or (e) a constant for a band this case needs (the Cys-hexose Amadori steps,
1-DG's exits, MGO's loss law)? Screening used titles and abstracts only (blinding holds).

**Tool and seeds.** `scripts/tools/openalex.py sweep` from `origin/openalex-tool` (PR #70), run from the scratchpad
through the wrapper, ~1 request per 1.1 s, anonymous, no e-mail.
- Forward and backward: Hofmann 1998 (10.1021/jf9705983; D1, D2, D4), Du 2023 (10.1016/j.foodres.2023.113337; D3),
  Zhao 2019b (10.1016/j.foodchem.2018.07.106; D3), de Roos 2005 (10.1021/bk-2005-0905.ch009; D5), Zhang 2026
  (10.1016/j.foodchem.2026.148681; part G constants), and, once the held-dossier check found them decisive-adjacent,
  Tressl 1994 (10.1021/bk-1994-0564.ch018; the positional label) and Voigt 2009 (10.1021/jf900459x; 1-DG's exits).
- Forward only: Meynier 1995 (10.1016/0308-8146(95)93282-v; D6), Xia 2022 (10.1021/acs.jafc.2c06639; D7), Kanzler
  2026 (10.1016/j.foodchem.2026.151175; no citers yet), Davidek 2006a (10.1021/jf060667q; 1-DG), Deng 2021
  (10.1021/acs.jafc.1c02023; the hexose ARP decay).

**Result.** 634 distinct works; 123 already in the registry by DOI; **511 unseen**, 315 with an abstract, 196 title only.
Every unseen work was assigned a decision (scratchpad table, not committed):

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
| not M4: a Vauthey sibling (Fanun, Leser, Aserin & Garti 2001, 10.1016/s0927-7757(01)00786-5), passed to H0's question | 1 |
| **candidates** (below) | 13 |

Notable exclusions, with the question they answer: Yaylayan & Keyhani 2000, 2001 (positional [13C]glucose with Ala /
Met, but pyrolysis GC-MS, no thiol at process temperature); Goto 1961 / Hayami 1961 (14C acetol, alkaline sugar
chemistry); Shibamoto & Russell 1976, 1977 (glucose-H2S-NH3, qualitative volatiles); Lee 2010 (GSH + glucose CAMOLA on
two thiophene aldehydes, not MFT; FD factors); Song 2019 ([13C6]glucose + Cys + 2-nonenal, MFT not named; M5); Yang 2015
(Cys + glucose with fat, M5); Cerny 2008, Schieberle 2005, Mottram 1994 (reviews of work already held).

**Candidates (each unseen; role fixed now, before reading; P = chance it holds the named measurement).**

| candidate (authorship from the OpenAlex record) | DOI | role, fixed now | the measurement that would count | P |
|---|---|---|---|---|
| Yeo & Shibamoto 1991, JAFC 39:370, microwave Maillard model system under different pH | 10.1021/jf00002a029 | **decisive for D2 (a second laboratory) only if** it is a glucose + Cys pot with MFT quantified at >= 3 pH; otherwise evidence | MFT at >= 3 pH, one design, stated scale | ~0.10 |
| El-Ghorab, Ashraf, Mohamed, El-Massry & Farouk 2010, Cys-ribose / Cys-glucose at pH 4.5 / 7 / 10 (no DOI in OpenAlex) | -- | as Yeo 1991 (D2, second laboratory) | as above; peak areas at fixed design suffice for a ratio | ~0.10 |
| Cannon & Ho 2023, J. Food Bioactives, key odorants in Arg / Cys-glucose | 10.31665/jfb.2023.18338 | **decisive for D1 (second laboratory) only if** NF and MFT are both quantified in one Cys + glucose pot; else part A evidence (Arg vs Cys) | NF and MFT in one pot, one assay each with stated recovery | ~0.10 |
| Du, Wang, Ma, Li, Bai, Cui, Fan & Wang 2024, LWT, 2-alkenals and Cys-glucose initial stage (BTBU) | 10.1016/j.lwt.2024.115741 | band source for the Cys-hexose Amadori steps (its no-aldehyde control), by amendment; never decisive (BTBU's pots are D3's laboratory) | a Cys-glucose Amadori time course with charges and T stated | ~0.3 |
| Schieberle, Fischer & Hofmann 2003, Weurman proceedings, CAMOLA (no DOI) | -- | evidence for D3 only (DFA); never decisive for D1 / D2 (same laboratory and possibly the same pots) | a glucose + Cys CAMOLA on MFT | ~0.2 |
| Tressl & Rewicki 1999, "Heat generated flavors and precursors" (chapter) | 10.1007/978-1-4615-4693-1_26 | **decisive for the which-carbon gap only if** it reports MFT from [1-13C]- or [6-13C]-glucose | MFT isotopologues from a positional glucose label | ~0.05 |
| Glomb, Gobert & Voigt 2010, ACS Symp. Ser. 1042, dicarbonyls from glucose and maltose | 10.1021/bk-2010-1042.ch004 | band source for 1-DG's exits and for `r_odg_dpo`'s ceiling (labelled sugars, O2 on / off); never decisive (no thiol) | a 1-DG -> C5 product (or its absence with a floor) | ~0.3 |
| Hayashi & Namiki 1986, Agric. Biol. Chem., sugar fragmentation in early browning | 10.1080/00021369.1986.10867692 | band source for R3's fragment supply pH law (C2 / C3 from glucose + beta-Ala negligible in acid, rising with pH) | fragment yields at >= 2 pH | ~0.5 |
| Pfeifer, Haase & Kroh 2013, JAFC, reactivity of thermally treated alpha-dicarbonyls (TU Berlin, Kroh group) | 10.1021/jf302959k | band source for MGO's self-loss law: a second laboratory beside Kanzler (lab boundary to be confirmed: Haase also co-authors Bork 2022b) | MGO alone at >= 2 temperatures | ~0.4 |
| Kanzler, Wustrack & Rohn 2021, JAFC, MGO + Ala or Lys melanoidins, pH 5 | 10.1021/acs.jafc.1c04594 | band source (MGO + Lys conversion; same laboratory as Kanzler 2026, so it widens no laboratory count) | MGO conversion vs time or T | ~0.3 |
| Yu, Xu & Yu 2016, JSFA, pH and amino acids on MGO in glucose-amino acid systems, 120 C | 10.1002/jsfa.8158 | part A evidence (Lys / Gly / Arg / Pro at fixed glucose and pH); band source for MGO supply vs pH | MGO, glucose and acid vs pH for >= 2 amines | ~0.4 |
| Luo, Xu, Li & Ho 2022, J. Food Bioactives, Gly / diGly / triGly ARP formation in a competing system | 10.31665/jfb.2022.18332 | part A evidence: peptide vs amino-acid competition for one sugar inside one study (D7; Ho network, as Xia 2022) | ARP ratios at fixed sugar across pH and T | ~0.5 |
| Luo, Zhu, Peng, Cui, Huang, Xu & Ho 2023, JAFC, Gly / diGly / triGly-glucose ARPs | 10.1021/acs.jafc.3c06250 | band source for peptide-ARP decay vs T (100 / 130 / 160 C) | ARP decay at >= 2 T | ~0.4 |

The combined chance that at least one D1 / D2 / which-carbon candidate holds a decisive measurement is
**~0.3** (1 - 0.9 x 0.9 x 0.9 x 0.95 = 0.31). None sets a band in this case until it is read and entered by an amendment
registered before any M4 engine result exists. The queue entries are in the private
`data/articles/DOWNLOAD_QUEUE.md` (gitignored), under "2026-11-05: M4 citation sweep".

**Held dossiers the sweep surfaced (already on disk, not in the 2026-11-04 case):** Tressl 1994, Voigt 2009, Voigt
2010, Smuda 2010, Du 2025, Zhang & Ho 1991b, Martins 2010, Schirle-Keller 1992, Hofmann 1997, Madruga 1995, Liu 2021.
Each now sits in sec. 2 or in 2.8 with its reason. A full screen of every held dossier that names glucose with a thiol
(~70 files, as M1c did for TTCA) was not done; open question 7.

## 3. The rivals

| id | part G: glucose + Cys -> MFT | part E: furfural supply vs pH (deferred) | part A: amines |
|---|---|---|---|
| R0 | **current engine**, scored twice: pinned at B60's constants, and with its own fitted coordinates refitted on the same rows as every rival: ~0 MFT (refused); FFT through `r_glc_fur` | `acid` = [H+] | one amine, folded |
| R1 | **the claim**: Glc + Cys -> Cys-hexose Amadori -> 1-DG + Cys -> 1-DPO + C1 (C-6 lost) -> the pentose steps (NF and DDP shares as M2 sets them) | amine-catalysed 1,2-enolisation, flat or weakly rising, plus an acid limb | classes with own constants, competition declared, catalysis of dicarbonyl loss separate from consumption |
| R2 | the same supply, entering **below norfuraneol**: 1-DG + amino acid (Strecker / reduction) -> a 1,4-dideoxy C6 (Huber 1990, CF-099) -> a 1,4-dideoxy C5 (DDP) -> MFT | M1's C1: the supply carries de Roos Fig. 4's measured U-shape (flat pH 3 - 5, rising above 5.5) and the split's pH law; no separate E8 term | -- |
| R3 | **C2 + C3 recombination** (B47's opened `k_glc_ha`), **now carrying Hofmann Table 6's measured HA + MGO -> NF** (NF 0.13 / 0.78 / 3.18 mol % at pH 3 / 5 / 7) and MGO's banded loss law (sec. 5.2) | the supply law is right; the steepness sits downstream (furfural loss, the HS- branch) | one constant per alpha-amino group, whatever the carrier |
| R4 | **the entry alone** on the existing amine-free 1-DG supply (the red team's band-now fix) | -- | -- |

**Notes on the rivals.**
- R1 and R2 both keep glucose C-1 as MFT's methyl (both lose C-6). Tressl 1994's only documented hexose C5 loses C-1
  (from 3-DG) and has no methyl. A route through a C-1-lost C5 that is then reduced to a 2-methylfuran is not proposed
  by any source on disk and is not a rival here; the positional experiment (sec. 4.4, arm 1) would expose it.
- R3 corrected (2026-11-05): the 2026-11-04 case read R3 as "no NF from fragments". Hofmann Table 6 measures HA + MGO
  -> NF at 0.78 mol % (pH 5), more NF than furfural, so R3 predicts fragment NF whose pH law rises (as Table 6). The
  engine has no HA + MGO -> NF step today; R3 gets one at Table 6's measured yields (a fed-intermediate constant from
  a table never run under the engine and not scored here, used as a fixed measured constant, not a band).
- R3 stays in the list although D3 rejects it by a measured share: the comparison reports that rejection, and R3 is
  the only rival whose MGO law Kanzler moves materially.

## 4. The decisive data

### 4.1 What cannot decide here

- **Seen pots** (sec. 0) can only test feasibility (development grade): Hofmann Table 1 (D4), the ribose half of D1,
  Hofmann Tables 3, 4, 8, 9, 10, Zhang 2026, de Roos Fig. 4, Bolton, Meynier, Deng 2021, Davidek 2006a, Xia 2022.
- **D3's pots** (Du 2023, Zhao 2019b) are M2's TEST: their route pattern is read from the published tables; no engine
  number on them before M2's item 7.
- **Band sources** (sec. 5.2) never decide, and decisive units never set a band.

### 4.2 The decisive rows

Directions follow from each structure. "[dry]" marks arithmetic on printed values.

| # | measurement (record) | value basis; engine status | R1 | R2 | R3 | R4 | separates |
|---|---|---|---|---|---|---|---|
| D1 | **NF / MFT in the glucose pot** (Hofmann 1998 Table 5: 7.1), with the ribose pot's 2754 as in-sample context | within-pot ratio of two printed values, no tau; glucose cells **run under R0 only**, ribose cells seen (fit rows) | glucose's 1-DPO feeds NF and MFT as ribose's does, so glucose NF / MFT equals the DPO-derived share's ratio that M2 sets; [dry] the measured NF ratio glucose / ribose is 2.5e-4 and the MFT ratio 0.096, 390x apart, so R1 needs most ribose NF to come from pentose-only routes (amine-free `r_pent_thermal_dpo`, `r_pent_caramel`) | glucose MFT bypasses NF: NF / MFT set only by glucose NF from other routes | NF from HA + MGO (Table 6) over MFT from HA + MP (Table 8): both rise with pH; its ratio follows MGO's standing level, hence Kanzler's law | as R1 | R1 from R2 (scale-free; nearest the intermediate). HexoseC5Entry's own docstring names this falsifier |
| D2 | **glucose MFT pH shape** (Hofmann Table 2: 0.3 / 1.9 / 0.4 at pH 3 / 5 / 7; ribose 55.3 / 19.8 / 2.5) | within-study contrasts log10(MFT_pH / MFT_5), no tau; pH 3 / 7 **run under R0 only**, pH 5 seen | low at pH 3 (2,3-enolisation minor in acid, P3), falling at pH 7 as beta-cleavage takes 1-DG (Davidek 2006a) | as R1 below 5; above 5 depends on the Strecker step's pH law | rising to pH 7, like fed HA + MP (0.01 / 0.23 / 0.27 mol %, Table 8) and fragment supply (Hayashi 1986, candidate) | the amine-free 1-DG supply's pH law, which no measurement on disk fixes | R3 from R1 / R2 by sign on the 5 -> 7 limb; the hexose shape against the pentose shape on the 3 -> 5 limb |
| D3 | **mixed isotopomers** (Du 2023, Zhao 2019b) | label share, no tau, no omega on the share (sec. 5.1); **reserved** | intact C5 | intact C5 | ~50 % mixed in a 1:1 pool | intact C5 | **rejects R3 by a measured share** (data, not a fit); does not separate R1, R2, R4 |
| D4 | **ribose / glucose MFT 10.4 and FFT 4.3** (Table 1) and their temperature trend: > 870 at ~100 C (de Roos Table I), 10.4 at 145 C | endpoint ratios; **seen** | follows the entry's barrier; [dry] if the two labs' pots were comparable, the hexose route's apparent barrier exceeds the pentose route's by at least ~128 kJ/mol (they are not: 37 % water lactate vs 0.5 M phosphate) | as R1 | follows MGO's loss law too: with Kanzler's Ea the C3 pool falls x13 - 33 faster at 145 C | the probe gave 45 - 78 at 145 C | side evidence only (seen; endpoints E6, E11 and the H2S budget also move) |
| D5 | **the Cys-hexose Amadori pool** (de Roos Table II; Zhang 2026 Fig. 3) with no hexose MFT at ~100 C (de Roos Table I, < 0.05 mg/mol Cys) | relative areas (II); a censored ratio > 870 (I); Table I glucose arm **not run**, Zhang seen | the pool exists and its decay to 1-DG is the slow step at 100 C | as R1 | no Amadori step needed | no Cys-hexose Amadori compound: the measured pool is unexplained | R1 / R2 from R4. de Roos I / II are screen signatures (sec. 5.5), so development grade after it; Zhang 2026 is evidence only |
| D6 | **FFT vs furfural vs pH** (Meynier Table 3; Hofmann F25#8) | within-study contrasts; seen | (part E, deferred) | | | | not separable on FFT alone; deferred with part E |
| D7 | **amine class at fixed sugar** (Xia 2022 TriG-ARP x5 - 7; de Roos Table II Fru-Cys / Fru-Gly 1 -> 4 - 5 with dilution; Brands 2001 casein leaves glucose loss unchanged; Kanzler 2026: amines catalyse MGO loss without consumption) | within-study ratios; Xia seen, de Roos II not run, Kanzler not run | each class its own constant and dilution law; catalysis declared | -- | (A-R3) equal per alpha-amine: Xia's x5 - 7 at equal molarity and de Roos's dilution-dependent ratio contradict it | -- | A-R1 from A-R3, within-study. One laboratory each |

**The plan's decisive data, checked.**
- "Hofmann ribose / glucose ratios (MFT 10.4, FFT 4.3)": D4 here, seen, side evidence.
- "The FFT pH slope": D6, deferred with part E.

### 4.3 The decisive data on disk, or their absence

- **R3:** rejected by D3's measured share; D2's 5 -> 7 limb agrees in sign. Not in question.
- **R4 against R1 / R2:** D5 (the pool exists) is the signature; it is screened (sec. 5.5), so the verdict there is
  development grade.
- **R1 against R2 (the claim against the bypass): one decisive unit, D1**, unseen under every rival but R0 (one
  laboratory, one temperature, NF by an unstated method, SIDA MFT). R0's existing numbers there are ~10 - 25 dex off
  (sec. 0) and cannot inform R1's or R2's bands, which come only from sec. 5.2. Whether it decides depends on how far apart R1 and R2 put glucose NF / MFT: under
  omega it needs >= 1.16 dex (x14) between the rivals' predictions (sec. 5.4). D2 does not separate R1 from R2 on disk
  (both predict the same sign on both limbs).
- **Which glucose carbon MFT loses, and any fed hexose intermediate -> MFT: absent.** No positional label on MFT, no
  fed 1-DG + H2S, no fed Cys-hexose Amadori compound -> MFT exists on disk. Four sweep candidates could supply D1 / D2 a
  second laboratory or the label (combined P ~0.3, sec. 2.9).
- **Part A:** one laboratory per class contrast; no Maillard flavour pot varies protein-bound vs free amine at fixed
  sugar, pH and temperature. Development grade at best.

**Therefore** R1 vs R2 is decided on disk only if the power check shows D1 separating them; otherwise M4's part G
verdict is development-grade ("R1 or R2, not R0, R3 or R4") and the decider is the designed experiment below.

### 4.4 The designed experiment (W6)

One laboratory, three arms. Throughout: n = 3 independent heatings per cell; sealed vessels with stated fill and
headspace (N2 and air arms in arm 2); pH stated on the cooled scale for every sample and in situ where an electrode
allows; every assay's LOD / LOQ and recovery reported; SIDA with labelled standards for MFT, FFT and NF.

| arm | conditions | species | decides |
|---|---|---|---|
| **1. which carbon** | Hofmann's design: Cys 33 + glucose 100 mM, 0.5 M phosphate pH 3 / 5 / 7, 145 C, 20 min; glucose as [1-13C], [6-13C] and 1:1 [13C6] / [12C6] (and unlabelled) | MFT, NF, FFT, furfural, MP, each isotopologue-resolved | MFT's methyl from C-1 with C-6 lost (R1, R2) vs C-1 lost (an unlisted 3-DG-side route); NF's label pattern equal to MFT's (R1) or independent of it (R2); mixed (R3). Also D2 and D1 in a second laboratory |
| **2. fed 1-DG** | synthetic 1-DG (Glomb 2000 route) 2 mM + Na2S or Cys 6 mM, +- glycine 33 mM, 0.5 M phosphate pH 5 and 7, 145 C, 2 - 20 min | MFT, NF, DDP and 1-DPO (OPD quinoxalines), acetic and formic acid, residual 1-DG | R1 (MFT without the amino acid, through 1-DPO and NF) vs R2 (MFT needs the amino acid; 1,4-dideoxy quinoxaline; no NF) vs R4 (amine-free 1-DG only) |
| **3. fed Cys-hexose Amadori** | synthetic D-fructose-L-cysteine (de Roos's preparation) 10 mM, 0.5 M phosphate pH 3 / 5 / 7, 100 and 145 C, 0 - 60 min | the Amadori compound, free Cys, 1-DG and 3-DG (OPD), MFT, NF | part G's supply constants (formation, decay, Cys return) and the mechanism of D2's pH maximum |

Request through W6 beside the BTBU / CAAS request for MFT per glucose consumed (ENGINE_VS_LITERATURE sec. 3.1, E1).

## 5. Tolerances, bands, screen and power

### 5.1 The measurement model per unit

| unit | role | value basis | reported noise | error terms in the primary | floor | density |
|---|---|---|---|---|---|---|
| D1 glucose NF / MFT (Hofmann T5) | decisive (R1 vs R2) | within-pot ratio | triplicates "differ by <= 10 %" (SIDA, MFT); NF method not stated: the same 10 % as a declared transfer, sensitivity x2 | omega per printed value (2 omega^2 on the ratio); no tau | MFT, NF: none printed; declared transfer 1 % of the lowest printed value in the table, required agreement over x0.1 - x10 (sec. 5.3) | exact log-normal; censored form exact |
| D1 ribose NF / MFT | context (in-sample) | within-pot ratio | as above | as above | as above | reported, not scored |
| D2 glucose MFT pH 3 / 7 against pH 5 (Hofmann T2) | decisive (R3 by sign; the pentose-shaped null) | two within-study contrasts sharing the pH 5 point (bivariate, correlation 0.5) | <= 10 % | omega per printed value; no tau | as D1 | exact bivariate normal in log10; censored where a rival predicts below the floor |
| D3 isotopomer shares (Du 2023, Zhao 2019b) | rejects R3 | label share | Du: no SD, two blind readers identical; Zhao: triplicate means | binomial-type share error from the printed precision (+-1 % absolute) plus a mixed-channel floor of 1 %; omega is a level error and does not apply to a share | M+1 - M+4 channels: 1 % (the printed "< 1" cells) | exact |
| D5 de Roos Table I (glucose < 0.05, ribose 43.5) | screen signature | censored ratio | not stated | omega per printed value on ribose; censored glucose | **0.05 mg/mol Cys (the paper's reporting limit)** | exact censored |
| D5 de Roos Table II | screen signature | relative RI areas (ratio) | response factors not given (mol ratios within ~x1.5, dossier) | omega per value plus the x1.5 response term | -- | exact |
| D5 Zhang 2026 MFT onset | evidence only (seen; a band source, so not screened) | ratios (RF 1) | n = 3, SD | omega per value | Zhang's printed LOD if present; else a declared transfer (1 % of the lowest printed value, required agreement) | exact censored |
| D7 class ratios (Xia, de Roos, Brands, Kanzler) | part A evidence | within-study ratios; Kanzler C/C0 | as printed; Brands digitised +-3 % of span; Kanzler +-1.5 % read, SD 0.5 - 10 % | omega per value; no tau | Kanzler: MGO <= 1 % of C0 (lowest calibrant) | exact |
| R0's sub-floor predictions | all units | a predicted level against a floor is an **absolute** comparison | -- | tau_ABS 0.964 dex on the predicted level, since it is an engine level against an assay floor | the unit's floor | exact censored |
| across laboratories (D4's 100 vs 145 C trend; Bolton's bound) | side | ratios in each lab | -- | omega per value; between-lab ratio offsets are not tau (tau is absolute-only) | -- | -- |

### 5.2 Bands and their sources

No band below comes from Hofmann Tables 2 or 5 (decisive), Du 2023 or Zhao 2019b (reserved), or from the screen's pots
(de Roos Tables I and II, Hofmann Table 1, Bolton). Starts are drawn inside the bands, never set from a decisive or
screened value. Some band sources have engine outputs in other records (Deng 2021, Martins 2003b, Davidek 2006a, Zhang
2026; sec. 0); none is scored or screened here. As in M1c, a seen source sets a centre or a fixed
measured constant, and any width it does not itself print is unsourced and carries required agreement.

| coordinate | band | source | kind |
|---|---|---|---|
| Glc + Cys -> Cys-hexose Amadori, formation | the trunk's Glc + Gly Schiff / Amadori constant (Martins 2003b class) x 0.3 - 3; Ea 101 kJ/mol +- 15 | Zhang 2026's Ea (BTBU, units not printed: the barrier only); the multiplier width has no source | same-class transfer; the width is **unsourced: required agreement** at both ends. de Roos Table II is a screen signature, so it does not set this |
| Cys-hexose Amadori -> 1-DG + Cys, decay | Deng 2021's Met-glucose ARP decay 0.0048 /min at 100 C (label pH 5.5), pH slope from its 7.5 arm, x0.3 - 3; Ea 85.7 - 107 kJ/mol | Deng 2021; Zhang 2026 (85.7); Martins 2005 (97 - 107) | same-class (hexose amino-acid ARP); the width unsourced: required agreement |
| Cys returned per Amadori decay | 0.71 - 1.0 | Martins 2003b, Davidek 2002 / 2005 (fed DFG: amine returned 71 - 96 %); Deng 2021 (~1) | measured, same class |
| `r_odg_dpo` (1-DG -> 1-DPO + C1) | 0 to a ceiling: the share of fed 1-DG not accounted for by the measured exits (Davidek 2006a; Voigt 2009; Smuda 2010) at their own conditions; barrier 40 - 135 kJ/mol | measured exits bound it from above only; Tressl reports no 1-DG C5 (no LOD) | numerical above zero: **a verdict that needs a value off the floor in this band does not ship without arm 2** |
| R2's Strecker / reduction step on 1-DG | 0 - the trunk's dicarbonyl + amino-acid Strecker constant | M3's class (once M3 has a verdict); Nedvidek 1992 shows the product exists | class bound; required agreement |
| NF and DDP shares downstream of 1-DPO | as M2 leaves them; if M2 has no verdict, M2's rivals carried as a factor | M2 | inherited |
| **MGO loss law (`k_mgo_mel` and any MGO sink on the network)** | Ea **95 - 121 kJ/mol**, level at 100 C pinned to B60's 0.0228 /min; order 1 (Kanzler's apparent order ~2 is a sensitivity) | Kanzler 2026 (80 -> 100 C, amino-acid pots) and Bork 2022 (MGO alone, 130 C; same laboratory) | measured class, **one laboratory, first reader only**: enters after its second reader; until then, and as a required-agreement sensitivity throughout, B60's Ea 20 |
| amine catalysis of MGO loss (part A) | excess k proportional to [amine]^0.75 - 0.84 (Ala); Lys x5 at 100 mM | Kanzler 2026 | one laboratory; part A evidence; no class constant fitted on flavour rows |
| R3's HA + MGO -> NF | Hofmann Table 6 yields at pH 3 / 5 / 7 | Hofmann 1998 Table 6 (not run, not scored; the same paper as D1 / D2 but other pots), used as a fixed measured constant | measured constant |
| R3's fragment supply pH law (`k_glc_ha` vs pH) | 0 - 1 dex per pH unit, rising | Hayashi 1986 (candidate, by amendment); until then numerical | numerical; required agreement |
| omega | 0.419 dex; sensitivity 0.390 | H0 (EX-059) | measured |

### 5.3 Floors

Hofmann 1998 prints no LOD for its SIDA thiols or its NF assay. The floor is a declared transfer (1 % of the lowest
printed value in the same table), with required agreement: every D1 / D2 verdict must hold for floors from x0.1 to x10 of
that value. de Roos Table I's "< 0.05 mg/mol Cys" is the paper's own reporting limit and is used as printed. Kanzler's
lowest calibrant (1 % of C0) is its floor. Zhang & Ho 1991b's MFT non-detect carries no LOD and is evidence, not a unit.

### 5.4 Power check (dry part done here; the synthetic part before the freeze)

**Dry arithmetic on the printed values, under omega 0.419 dex per value** (one-degree-of-freedom threshold 1.92 nats,
which a favoured rival exactly at the datum reaches only if the other rival misses by >= 1.96 sd):
- **D1 (glucose NF / MFT, two values):** sd = 0.419 x sqrt(2) = 0.593 dex. It decides between two rivals only if their
  predicted glucose NF / MFT differ by >= 1.16 dex (x14). If one rival predicted equal NF / MFT ratios in the glucose
  and ribose pots (the naive R1), the ratio-of-ratios form (four values, sd 0.838 dex) would put the datum (-2.59 dex)
  3.1 sd away, ~4.8 nats; R1 does not predict that, because the ribose NF has pentose-only sources.
- **D2, 5 -> 7 limb:** measured -0.68 dex; R3's transfer from Table 8 +0.07 dex; difference 0.74 dex = 1.25 sd, **~0.8
  nats: below the threshold.** D2 alone cannot reject R3; D3 does.
- **D2, 3 -> 5 limb:** measured +0.80 dex against the ribose (pentose) shape -0.45 dex: difference 1.25 dex = 2.1 sd,
  ~2.2 nats: it rejects "the hexose route has the pentose route's pH law" marginally (one limb; the shared pH 5 point
  makes the bivariate form a little weaker).
- **D3:** a 1 % mixed channel against R3's ~50 %: decisive at any omega (a share, not a level).

**Synthetic part, required before the freeze.** On the first build of the network rivals (W3), with every rival at
its band centre and over its bands (Monte Carlo, the sec. 5.1 densities), simulate D1 and D2 from each rival, score all
rivals, and report for each pair P(the generating rival wins by > 1.92 nats). A unit whose power is < 0.5 for a pair is
demoted to "ranks, does not decide" for that pair before any measured decisive value is compared with a rival. The check
uses synthetic data only; R0 at its shipped constants may be evaluated on the designs (as the M2 prereg did).

### 5.5 The mechanism-signature screen (before any fit; demotions only)

Plug-in runs at band centres and corners, on the one network, on signatures from pots that are already seen or that
this case gives up as decisive:

| signature | pot | passes if |
|---|---|---|
| S1 the hexose route shut at ~100 C | de Roos Table I (glucose < 0.05 vs ribose 43.5 mg/mol Cys; lactate pH 5, label) | predicted ribose / glucose MFT > 870 (censored) |
| S2 a Cys-hexose Amadori pool forms in situ | de Roos Table II (Fru-Cys / Fru-Gly 0.9 - 1.1 at 15 % water) | the pool's predicted ratio within the x1.5 response term and omega |
| S3 consistency with the shipped T1 ratio | Hofmann Table 1 (ribose / glucose MFT 10.4, pH 5) | within omega on the ratio |
| S4 Bolton's bound | glucose bouillon, 120 C | predicted Cys-derived MFT <= the bound |

Zhang 2026's MFT onset (absent before 90 min) is **not** a screen signature: its Table 2 barriers set bands (sec. 5.2),
and one paper's pots cannot both set a band and screen it.

After the screen, a rival that fails a signature inside all its bands is demoted (reported, not dropped); none is
added, promoted or redesigned around the screen. The screen's pots are then declared seen; D1 and D2 are untouched.
Stop rule: if no rival among R1, R2, R4 passes S1 and S2 inside its bands within ~15 h, part G reads "hexose entry not
identifiable on disk" and the designed experiment goes out.

## 6. What M4 needs from the one network (W3)

- **Species:** Glc, Fru, Cys, the Cys-hexose Amadori compound (new), 1-DG shared by the hexose and sulfur modules,
  1-DPO, DDP, NF, MFT, furfural, FFT, HA, MGO, MP; for R2 a 1,4-dideoxy C6. Today 1-DG lives in the trunk and 1-DPO in
  the sulfur lane, which do not compose.
- **One MGO pool** with every sink declared (`r_mgo_mel` with a banded temperature law, `r_mgo_mp`, `r_cys_actz`,
  `r_mgo_ak`), and an HA + MGO -> NF step for R3 (Hofmann Table 6).
- **Amine classes as species, not a convention:** free alpha-amino acid (with identity factors), cysteine,
  peptide (GSH, diglycine), protein-bound lysine (one key, A2b), each declaring which sugars it competes for, so
  `compose` refuses undeclared competition; and a separate catalytic term (amine catalysis of dicarbonyl loss without
  consumption, Kanzler 2026).
- **One pot with all of them:** Glc + Cys + Gly + protein in one integration (today the lanes refuse to combine).
- **In-situ pH on the hexose Amadori and enolisation steps,** with the Cys-hexose Amadori compound's pKa in the one
  titratable table; the cooled-scale observable for label-pH comparisons; and a per-step hook for M5's buffer factor
  (every Hofmann pot is 0.5 M phosphate).
- **Isotopologue pools:** labelled and unlabelled copies of the sugar carbon, position-resolved at C-1 and C-6 (or atom
  tracking), so D3 and arm 1 are observables, not hand calculations.
- **The forced-input driver** (A2 contract sec. 2): 1-DG / 3-DG courses (Deng 2021) clamped while the entry is tested.

## Open questions for the build

1. **Order with M1 and M2.** Part G feeds the pentose steps whose NF / DDP shares M2 refits and whose supply M1's C1
   replaces. D1 is read on the network after M2's verdict, or with M2's rivals carried as a factor.
2. **The Hofmann glucose pH 3 / 7 arms.** DEV under both frozen splits, but declared Maillard-path hold-out bundles.
   This case treats them as decisive (plan sec. 7: held-out status is the split files'), scored once in item 6, never
   a fit row. If the older declaration binds, D2 moves to item 7 and R1 vs R2 rests on D1 alone.
3. **The Cys-hexose Amadori step's constants.** Zhang 2026's Ea ~101 is from a GSH pot with units not printed; de Roos
   Table II gives only relative areas and is a screen signature. The formation multiplier's width is unsourced (required
   agreement); Du 2024b (candidate) may source it by amendment.
4. **E8 after M1.** If C1 ships, does part E survive for `r_glc_fur` and `r_arp_tdp` alone? Decide before building
   the amine-catalysed term. Schirle-Keller's 147 kJ/mol furfural law (glucose + Cys, unbuffered) bears on `r_glc_fur`.
5. **Part A's constants** come from one laboratory per class (Xia, de Roos, Brands, Tang, Kanzler). Declare them as
   bands with a per-class basis; no class constant is fitted on flavour rows.
6. **Kanzler's second reader** and the CF entry for `k_mgo_mel` (proposed text in the dossier, sec. 6 item 1). Until the
   second reader, the MGO band is a sensitivity only, and R0 is scored at Ea 20 regardless.
7. **A full held-dossier screen** (every dossier that names glucose with a thiol or 1-DG, ~70 files; two blinded
   readers, as M1c did) before the freeze? The sweep found eleven held sources the first draft missed, Tressl 1994
   the most important.
8. **W6 request:** the BTBU / CAAS MFT yields per glucose consumed (ENGINE_VS_LITERATURE sec. 3.1, E1) would give D4
   a second laboratory; the designed experiment (sec. 4.4) would decide R1 vs R2.
9. **D1 and D2 after R0's probe.** Both pots were run under R0 (the no-refit thiol sub-model probe) before this case,
   and the author read that probe's per-data-set summary while refreshing. This case keeps them decisive between R1,
   R2, R3 and R4, which have never run there. If the review rules that any engine number on a pot spends it, part G has
   no unseen decisive unit on disk and its verdict is development grade until the designed experiment or a sweep
   candidate (Yeo 1991, El-Ghorab 2010, Cannon 2023) supplies one.

## Revision log

- **2026-11-04 (v1).** Items 1 - 5 written before any result; no engine run.
- **2026-11-05 (refresh), not frozen.** Still before any result; no engine run on any pot; the sweep read titles and
  abstracts only.
  - This week's rules applied: value basis per unit, tau on absolute levels only, omega 0.419, floors from assays,
    exact densities, the adequacy rule, R0 pinned at B60, bands and starts from unscored sources with required
    agreement for unsourced widths, the screen before any fit (demotions only), the power check, pH scales named.
  - Sec. 0 declares which pots have been seen under the engine; Hofmann Table 9 and Zhang 2026 moved to seen; D1's
    and D2's glucose pots are declared run under R0 only (the no-refit thiol sub-model probe, whose per-data-set
    summary the author read during the refresh).
  - Kanzler 2026 added (fed MGO; the MGO sink's Ea 95 - 121 vs `k_mgo_mel`'s 20 at its floor) to sec. 1, 2.2, 2.5, 3,
    5.2 and 6.
  - Held sources added: Tressl 1994 (positional labels: the "no positional label on disk" line was wrong; none is on
    MFT), Voigt 2009 / 2010, Smuda 2010, Du 2025, Zhang & Ho 1991b, Martins 2010, Schirle-Keller 1992; Hofmann 1997,
    Madruga 1995, Liu 2021 left out with reasons.
  - Hofmann Table 8 corrected (HA + MP, not C2 + C3 + H2S); Table 6 (HA + MGO -> NF) added, and R3's D1 prediction
    corrected from "no NF from fragments".
  - The citation sweep (sec. 2.9): 12 seeds, 634 works, 511 unseen, 13 candidates with roles fixed before reading.
  - D1's primary is now the glucose pot's NF / MFT (unseen), with the ribose ratio as in-sample context.
  - The designed experiment named (sec. 4.4): positional labels with MFT and NF, fed 1-DG, fed Cys-hexose Amadori.
