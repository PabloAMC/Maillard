# Misfit inventory: every open structural misfit between the engine and the data

> **The misfit protocol's inventory** ([OPERATING_PROCEDURE.md](../guides/OPERATING_PROCEDURE.md) sec. 4). Hand-kept, like the
> registers. Compiled 2026-11-06 from `origin/main` at `72ce2db4`, plus two unmerged branches: `m3-phase1-build`
> (M3 phase-1 records) and `decisive-measurements` (PR #79). Read-only: no fit, no engine change, no new score. No
> register, case, prereg, `src/` or `results/` file was edited to write this page.
>
> **What it is for.** The protocol has so far been applied one misfit at a time. This page puts every open structural
> misfit in one place, for three reasons: so misfits that share a cause are seen together; so each leading candidate
> mechanism gets its cheapest check on data already seen; and so those checks run in batches, where one check that
> tests a shared piece replaces several separate ones.
>
> **Spent data: the rule every case must respect.** A check *spends* every dataset it reads a value from. A dataset
> spent on a candidate cannot later be that candidate's decisive test, nor the decisive test of any rival that shares
> the piece. Each check below lists what it spends (column *spends*). A case that later names one of these datasets as
> decisive must list it as SEEN under its item 2.
>
> A check marked **NOT RUNNABLE** needs a blind or reserved dataset. Do not run it, and do not read that dataset to
> scope it. A check marked **owned by \<case\>** reads data that an open or frozen case has named as decisive. It runs
> inside that case's pre-registration or not at all.

Conventions:
- **dex** = log10(predicted / measured): + means the engine reads high. Whitened residuals are in σ, as their record
  prints them.
- **Protocol steps:** 1 integrity, 2 conditions, 3 measurement, 4 mechanism, 5 unresolved. "Step reached" is the
  highest step with recorded work. "Never checked" lists the lower steps with none.
- **Integrity has two levels.** **K** means the rows were keyed by two readers. **R** means the source was re-read
  *because of this misfit*. Only R counts as step 1 done.
- **W5.4 cells** are odorant × lever cells of the lever table ([STATUS_AND_PLAN.md](../guides/STATUS_AND_PLAN.md)
  W5.4).
  - Odorants: **Th** meaty thiols (MFT, FFT, H2S), **St** Strecker aldehydes, **Pz** pyrazines, **La** lipid
    aldehydes, plus **AA**, acrylamide, the safety output.
  - Levers: **cys** (cysteine or yeast extract), **sug** (sugar type), **thi** (thiamine), **pH**, **fat**, **Tt**
    (temperature and time), **ck** (cooking method).
- **Decisive measurement** links to `docs/reference/decisive_measurements.yml` on PR #79 (branch
  `decisive-measurements`, not yet on main) as `DM:<id>`. Where no entry exists, it names the CF decider or the atlas
  design ([ENGINE_VS_LITERATURE.md](pathway_atlas/ENGINE_VS_LITERATURE.md) sec. 3, `E#`).
- **Pass rule (every check; draft).** A check **EXPLAINS** if, with constants from cited measurements only, the prediction
  falls within the unit's own noise band of the datum without moving in-sample rows out of theirs. The table states
  each check's *unit*, *noise band* and *in-sample rows*. The in-sample rows are the B60 objective rows of the touched
  lane, unless the check names others. A check that needs a fitted constant to pass does not explain; it records
  "sufficient only if \<constant\>". A check **REFUTES** a candidate when the candidate's best cited-constant
  prediction misses the datum by more than twice its noise band on the decisive side.

---

## 1. Summary

**42 open structural misfits.** The counts below are by lane and by highest protocol step reached. In the sulfur lane,
"step 4" mostly means mechanism candidates were tested and none was separated (EX-061, EX-062, cluster runs).

| lane | n | step 1 | step 2 | step 3 | step 4 | step 5 | step 1 (R) never done |
|---|---|---|---|---|---|---|---|
| sulfur (thiols, furfural, NF, TTCA, pH balance) | 18 | 1 | 1 | 1 | 14 | 1 | 4 |
| trunk (dicarbonyls, acids, glycation, buffer) | 9 | 2 | 0 | 0 | 6 | 1 | 2 |
| Strecker (engine B18 / B22 / B29; M3 layer) | 5 | 2 | 1 | 1 | 1 | 0 | 4 |
| pyrazine | 3 | 0 | 0 | 0 | 3 | 0 | 1 |
| acrylamide (incl. furanics) | 6 | 4 | 0 | 0 | 1 | 1 | 4 |
| lipid | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| **total** | **42** | 10 | 2 | 2 | 25 | 3 | **16** |

How the steps are assigned:
- "Step 5" is a misfit that is recorded and refused by name, with its alternatives carried: MF-14 (H0), MF-40 (the
  low-T acrylamide guard) and MF-20 (the decisive measurement only).
- A misfit with no recorded work at any step is counted under step 1, with a gap.

---

## 2. The inventory

Entry fields:
- **Obs**: the observable and its lane.
- **Size**: the size, with the record it comes from.
- **Sig**: the structural signature.
- **Step**: the step reached, with its evidence, and what was never checked.
- **Cand**: the candidates on record.
- **Shared**: the shared structure.
- **Checks**: section 3.
- **W5.4**: the lever-table cells affected.
- **Decisive**: the decisive measurement.

### Sulfur lane

**MF-01 — The MFT / FFT time shape: no late rise, peak too early**

- **Obs:** MFT ratio rows in ribose + Cys, 0.5 M phosphate pH 5, 100 C (Schieberle 2000 Table IV); the Nestle
  110 C curve; the TTCA ladder T-order.
- **Size:**
  - MFT 360/30 −1.09 dex, 720/30 −1.42 dex at B60 (was −1.80 at B9). The R0 peak is at 158 min, while the pot is still
    rising at 12 h. Source: `kinetic_core_m1_preflight_checks.md` D1.
  - The R3En flag is −1.04σ over 8 rows.
  - B60 gives the Nestle curve 0.42 / 0.80 / 1.05 / 1 against 0.008 / 0.067 / 0.31 / 1 measured (CF-015).
  - Zhai 2023 JAFC thiols order 100 < 120 < 140 C in 16 / 16 comparisons; the engine gets it right in 1 / 8 (CF-017).
- **Sig:**
  - right level, wrong shape;
  - one sign from B9 to B60;
  - `k_mft_decay` frozen, a constant no source pins (CF-107);
  - `Ea_decay_thiol_sink` on its 102 kJ/mol ceiling.
- **Step 4:**
  - Integrity R: CF-017 (re-read, error confirmed); CF-055 (Schieberle text vs table).
  - Conditions: the oxidant is ruled out (preflight G1: OX use ≤ 0.009 %).
  - Measurement: the M1 review found three major measurement-model faults (`kinetic_core_m1_results_review.md`).
  - Mechanism: EX-061 R0 vs C1 is **not separated**; the two miss by about 1 dex in opposite directions.
  - Never checked: none.
- **Cand:**
  - (a) the sink is too fast: US 10,076,126's co-heated labelled standards survive, so the sink is ≥ 5x too fast
    (CF-013);
  - (b) slow supply (C1, M1);
  - (c) disulfide reservoir release (C3; Guth 1995).
- **Shared:** `ch_mft_decay` / `k_mft_decay`, `Ea_decay_thiol_sink`, `k_dimer_release`, supply (MF-02).
- **Checks:** CK-23, CK-24, CK-25, CK-27.
- **W5.4:** Th × Tt, Th × ck, Th × cys.
- **Decisive:** DM:thiol_supply_structure; atlas E12 (labelled MFT spiked mid-cook).

**MF-02 — Cysteine and pentose consumed in seconds**

- **Obs:** free Cys and pentose courses (Zhang 2025b blank, Cao 2017, Zhai 2019).
- **Size:**
  - Zhang 2025b Cys half-life: 5.13 min measured, 0.52 s in R0 (≈ 2.8 dex in time).
  - Cao pentose at 300 min: 1.09 measured vs 0 in R0.
  - C1 overcorrects: Cys −0.71 dex, PENT −0.52 dex in rate.
  - Source: preflight D2; M1 review.
- **Sig:** one sign across 3 laboratories (BTBU, DFA, Ho network); wrong shape.
- **Step 4:**
  - Integrity R: CF-063 (Zhang's single ELSD curve) and CF-109 (points below the curve's zero), so supply credit rests
    on Zhang alone.
  - Mechanism: M1 (EX-061) and the supply screen (EX-062: no candidate feasible on all three signatures).
- **Cand:**
  - (a) a TTCA / CAMA reservoir (`TTCA_FORMATION`, `CYS_AMADORI_ROUTE`);
  - (b) Zhang's low points are censored (CF-109).
- **Shared:** `r_ttca_cys`, `r_cama_form`, `r_cys_thermal`, `r_pent_caramel`; G4.
- **Checks:** CK-30 (owned by M1c).
- **W5.4:** Th × cys, Th × sug, Th × Tt.
- **Decisive:** DM:thiol_supply_structure.

**MF-03 — The TTCA exit has no pH law; the TTCA quotient's pH dependence is unmatched**

- **Obs:**
  - fed-TTCA loss, Zhai 2021 Fig. 4b (pH 5.5 → 8);
  - K_app = [TTCA] / ([Cys][Xyl]), Zhao 2017b vs Cao 2017;
  - Zhang 2025b's plateau.
- **Size:**
  - Loss ratio pH 8 / 5.5: 9.60 measured vs ≤ 1.36 for every candidate (≥ 0.85 dex short).
  - Zhang's plateau needs K_app ≥ ~19 L/mmol, against Cao's 0.060 L/mmol (≈ 2.5 dex).
  - Zhao's slope is −0.29 dex/pH against Cao's flat +0.01.
  - Sources: `kinetic_core_supply_structure_signature_screen.md`; CF-110; CF-111.
- **Sig:** correlates with pH; leave-one-group-out tension (Cao vs Zhang vs Zhai).
- **Step 4:**
  - Integrity R: CF-010 (pH of Zhai 2021, fixed) and CF-108 (printed Ea unusable).
  - Conditions: in-situ pH makes the needed exit steeper, not shallower (CF-110 expl. 2).
  - Mechanism: screen signature (b) fails every candidate.
- **Cand:**
  - (a) a base-favoured exit (deprotonated ring N, or a base-catalysed Amadori step);
  - (b) K(pH) through Kallen's speciation;
  - (c) read error in Zhao's late low-Cys points.
- **Shared:** `REACTION_PH_FACTOR` (no entry for `r_ttca_deg`, `r_cama_form`, `r_ttca_cys`); the in-situ pH state
  (G1).
- **Checks:** CK-31 (owned by M1c), CK-12 (side reading).
- **W5.4:** Th × pH, Th × cys.
- **Decisive:** DM:thiol_supply_structure, arms 1 and 3.

**MF-04 — TTCA exchange release 36–140x too fast**

- **Obs:** cysteine released from TTCA in the dilute fed pots (Zhai 2021) vs the concentrated pots (Zhai 2019, Cao).
- **Size:**
  - `k_ttca_cys` is 0.0286 /min at 100 C, against the feasible ridge of 0.9–3.2e-4 (≈ 2 dex).
  - The engine releases Cys 2–30x faster than TTCA is lost.
  - Source: CF-076.
- **Sig:** constants on bounds (three bounds active at the feasible corner); leave-one-group-out tension, inside one
  laboratory.
- **Step 4:**
  - Integrity R: CF-091 (Fig. 7a is a difference of two pots) and the kallen1971 dossier correction.
  - Mechanism: cluster 5 and resolve2_sulfur. Exchange E0 is REJECTED in two laboratories. E4c is admitted
    provisionally (P 0.45), feasible only on a (k25, Ea) ridge.
- **Cand:**
  - (a) the in-situ-pH Kallen law with a low release Ea (P 0.25);
  - (b) a slow exchange with a large pH-dependent K (P 0.2);
  - (c) a second Cys sink, or an unmeasured adduct pool (P 0.2).
- **Shared:** `r_ttca_cys`, K(pH), the TTCA ring-N pKa; G4, G1.
- **Checks:** CK-32 (owned by M1c; seen data do not separate the candidates).
- **W5.4:** Th × cys, Th × pH.
- **Decisive:** DM:thiol_supply_structure. The CF-076 decider is one Arrhenius series of Cys release at pH 4.5 and 6.

**MF-05 — Cys-Amadori formation: Xie phosphate pots 3–6x faster than Jiangnan water pots**

- **Obs:** CAMA formation rate per TTCA (Cao 2017, Zhang 2025b vs Zhai 2019).
- **Size:**
  - The joint-fit residual gap is 3–6x (0.5–0.8 dex).
  - Cao's Cys-Amadori rows are −0.33 to −0.57 dex.
  - K sits on its upper bound 0.075.
  - Source: CF-011.
- **Sig:** one sign by buffer class; constant on a bound.
- **Step 4:**
  - Integrity: Zhai's ELSD calibration for CAMA was never printed (an author request).
  - Conditions: monotone pH and concentration laws are excluded (cluster 5).
  - Mechanism: phosphate catalysis is sufficient (DEV κ_P [5, 40] /M) but not shown.
- **Cand:**
  - (a) phosphate catalysis of formation (P 0.5);
  - (b) Jiangnan ELSD quantification of low CAMA (P 0.3).
- **Shared:** a per-step buffer factor (G5); `r_cama_form`.
- **Checks:** CK-33 (owned by M1c / M5; no seen-data discriminator).
- **W5.4:** Th × ck (buffer and matrix), Th × cys.
- **Decisive:** CF-011 buffer arm C (water vs 0.05 / 0.2 M phosphate, plus citrate at fixed pH).

**MF-06 — The cooled-pH balance: too acidic in buffered Cys pots, wrong shape in unbuffered ARP pots**

- **Obs:** cooled pH after heating. Meynier 1995 (negative control), Hou 2026, Zhou 2023 Fig. 2, Cao 2017 Fig. 4c.
- **Size (pH units, not dex):**
  - Meynier: engine drift −0.79 / −0.65 / −0.39 against ≤ 0.2 measured.
  - Hou: mean +0.58, RMS 0.69 (tolerance 0.3).
  - Zhou label sensitivity is about 4x too large.
  - Cao is −0.20 at 30 min.
  - The acid yield sits ON ITS LOWER BOUND in B60.
  - Source: CF-084.
- **Sig:**
  - one sign across 3 laboratories;
  - constant on a bound;
  - calibrated on one arm (Zhou 2023 ARP + Cys, three anchors).
- **Step 4:**
  - Integrity R: CF-084 note (4) on the cooled scale.
  - Conditions: cluster 3.
  - Mechanism: outside support for candidate (a) from Sohn 1995 (NH3 leaves cysteine as a base, 76–82 %) and Hou
    2026 (glycine released, 0.25–0.40 per ARP lost).
  - Candidate (c) alone was tried (pKa 7.07 / 7.6 / 8.18, RMS 0.69 / 0.45 / 0.84): it fails.
- **Cand:**
  - (a) the amine is returned, not destroyed (`AMINE_FATE_BASIS`, `CENTRE_LEDGER`);
  - (b) acid forms with a lag downstream of the ARP;
  - (c) the ARP pKa.
- **Shared:** `ph_state` (G1). Every pH-dependent rate downstream of it.
- **Checks:** CK-12, CK-13, CK-14.
- **W5.4:** pH × every odorant (the pH lever is read through this balance in unbuffered foods).
- **Decisive:** the CF-084 remedy (a measured-trajectory operator). No DM entry.

**MF-07 — Norfuraneol decay is pH-free and too slow**

- **Obs:** fed-NF loss. Bork 2022 Fig. 1D and Bork 2023 Fig. 4B (pH 5, 130 C); Mikami 2017 / Nakamura 2020 (pH 8,
  100 C).
- **Size:**
  - At pH 5 / 130 C the engine is 3.1–4.3x too slow (−0.5 to −0.63 dex in k; CF-112).
  - At pH 8 it is ~700x too slow (F13-a).
- **Sig:**
  - one sign in two laboratories at two pH values;
  - correlates with pH;
  - `k_nf_decay` is unpinned (`unpinned_constants.yml`).
- **Step 4:**
  - Integrity R: CF-112 / CF-113, two readers, within-lab spread x1.4.
  - Mechanism: CF-117 (metal / O2 dependence at pH 8 in Nakamura's text, not quantified); the B61 hydroxide term is
    off.
- **Cand:**
  - (a) the pH-free constant is too slow at acid pH;
  - (b) a metal / O2-catalysed enolate route (P1');
  - (c) assay or vessel effects (the x1.4 spread).
- **Shared:** `r_nf_decay`, `NF_HYDROXIDE_DECAY`; G9 (O2).
- **Checks:** CK-39 (owned by M1b(ii)).
- **W5.4:** Th × pH, Th × Tt.
- **Decisive:** DM:norfuraneol_decay_ph_law.

**MF-08 — Engine norfuraneol at pH 3 is 55x below pH 5 (CF-012 residual)**

- **Obs:** NF in the Hofmann 1998 ribose + Cys pH ladder, 145 C, 0.5 M phosphate.
- **Size:**
  - Engine NF pH 5 / pH 3 is 55x (1.74 dex), i.e. ~1.4 dex per in-situ pH unit.
  - Against this, Meynier's NF *falls* with pH (wrong sign).
  - Source: CF-012 tests.
- **Sig:** correlates with pH; wrong sign against a second laboratory.
- **Step 4:**
  - Integrity R: CF-012.
  - Conditions: cluster 3 shows the in-situ contrast is pH 3.71 vs 5.27 → 4.47; the engine's pH-5 pot probably runs
    too acidic (CF-084).
- **Cand:**
  - (a) the NF formation acid law (`r_dpo_nf` via the pentose split);
  - (b) the in-situ pH (G1);
  - (c) the amine fate.
- **Shared:** G1, G2.
- **Checks:** CK-15.
- **W5.4:** Th × pH.
- **Decisive:** DM:pentose_split_vs_removal.

**MF-09 — The pentose → furfural / FFT supply pH law is too steep; MFT / FFT moves the wrong way**

- **Obs:** FFT, furfural, MFT and ACTZ pH ladders. Hofmann 1998 T2, Meynier 1995, Cerny 2007, Xu 2007.
- **Size:**
  - FFT pH 4.5 → 6.5: x400 engine vs x6.1 measured (+1.82 dex too steep).
  - FFT pH 3 → 7: x3600 vs x19 (+2.28 dex).
  - Xu 2007 log ratios, predicted vs measured: FUR +2.62 vs +0.26, MFT +2.35 vs +0.80.
  - MFT / FFT goes 0.40 → 3.1 in the engine against 1.79 → < 0.07 measured: opposite sign.
  - Source: ENGINE_VS_LITERATURE sec. 2.0, App. A; `kinetic_core_contrast_benchmark_score_v1.md`.
- **Sig:** one sign across 3 laboratories; correlates with pH.
- **Step 4:**
  - Mechanism: E8 is SUPPORTED; candidate laws are in the atlas.
  - Never checked: the attribution between supply and step inside the engine (CK-18).
- **Cand:**
  - (a) the `acid` = [H+] factor on `r_pent_tdp`, `r_arp_tdp` and `r_glc_fur` should be amine-catalysed, flat or
    rising with pH (hexose k4 +0.35 dex/pH, F06-b);
  - (b) the steepness sits in the FFT step's sulfide speciation, not the supply.
- **Shared:** G2 (the enolisation pH law); the HS⁻ fraction (G10).
- **Checks:** CK-17, CK-18.
- **W5.4:** Th × pH, Th × sug.
- **Decisive:** atlas E8 (ribose + Cys, pH-stat 4.5 / 5.5 / 6.5, 140 C).

**MF-10 — MFT runs ~100 % through free norfuraneol**

- **Obs:** route shares of MFT. Cerny 2003 pot D, Givaudan, Hofmann 145 C.
- **Size:**
  - NF share 0.999 vs 0.07 (+1.15 dex).
  - Givaudan ribose / NF MFT 0.0066 vs 48 (−3.86 dex).
  - Hofmann 21 vs 2.6 (+0.91 dex).
  - Source: atlas E6.
  - Pot D is 6.9 sd under ω in M2.
- **Sig:**
  - one sign across 3 laboratories;
  - `k_nf_mft` (−2.67) and `k_ddp_mft` (−6.54) 3.9 dex apart;
  - the R3En offset on `k_nf_mft` moves the wrong way.
- **Step 4:**
  - Integrity K.
  - Mechanism: the atlas red team (E6 strongly supported at ≤ 100 C; it needs E11).
  - Never checked: R for the Givaudan ratio.
- **Cand:**
  - (a) the direct DDP route is under-weighted (`r_ddp_mft`);
  - (b) one barrier makes the share T-invariant (E11, MF-11).
- **Shared:** `r_nf_mft`, `r_ddp_mft`, `r_dpo_nf`, `Ea_lumped_formation`; G7.
- **Checks:** CK-40 (owned by M2), CK-26 (side reading).
- **W5.4:** Th × sug, Th × Tt.
- **Decisive:** atlas E6 (Cerny pot D at 120 / 145 C).

**MF-11 — One formation barrier; the thiol-sink barrier sits on its ceiling**

- **Obs:** temperature dependence of the thiols.
- **Size:**
  - `Ea_lumped_formation` is 64.08 kJ/mol, against measured classes 81–147.
  - `Ea_decay_thiol_sink` is 101.25 against the Gigl ceiling of 102: the 4th consecutive wave on its ceiling (C1 also
    binds there).
  - CF-017's T-order: 1 / 8 right.
  - The stale lane offset tracks temperature (ρ −0.82; `lane_offset_diagnostic.md`, pre-H0).
- **Sig:** constant on a bound and far from its measured analogue; correlates with temperature.
- **Step 4:**
  - Integrity R: the CF-017 slice artefact is removed.
  - Mechanism: atlas E11.
  - Never checked: forward runs with the class barriers.
- **Cand:**
  - (a) per-route measured barriers;
  - (b) the sink Ea compensates instant supply (CF-107).
- **Shared:** `formation_Ea_by_route`, `FURFURAL_BARRIERS`, the thiol-sink decay family; G7, G3.
- **Checks:** CK-26, CK-27.
- **W5.4:** Th × Tt, Th × ck.
- **Decisive:** atlas E11 (one pentose + Cys pot at 100–160 C with route markers).

**MF-12 — Hexose gives almost no MFT**

- **Obs:** glucose + Cys MFT and FFT (Hofmann 1998 T1; hold-out bundles at pH 3 and 7).
- **Size:**
  - Ribose / glucose MFT: 4.4e10 vs 10.4 (−9.6 dex on glucose MFT).
  - FFT ratio: 0.47 vs 4.3 (−0.96 dex).
  - 8 panel rows are refused.
- **Sig:** `k_glc_ha` bound-limited (its admissible band is the whole declared band); wrong route (100 % C2 + C3 vs
  ≤ 6 % by labelling).
- **Step 4:**
  - Integrity R: atlas sec. 6 recount.
  - Mechanism: the red team REFUTED the band-now fix (`HEXOSE_C5_ENTRY` gives 45–78, or 0.79 with Gly). A new
    Glc + Cys → 1-DG step is needed.
- **Cand:**
  - (a) a Glc + Cys → 1-DG step;
  - (b) the hexose C5 entry plus E6.
- **Shared:** `r_glc_c2c3`, `r_odg_dpo`, `HEXOSE_C5_ENTRY`; G8.
- **Checks:** CK-41 (owned by M4; seen-data step done).
- **W5.4:** Th × sug.
- **Decisive:** DM:glucose_mft_carbon_origin.

**MF-13 — Fed furfural / NF + sulfide: the engine converts everything in seconds, then destroys it**

- **Obs:** FFT from furfural + H2S and MFT from NF + H2S (Hofmann 1998 T3 / T4; Munch 1997 T5).
- **Size:**
  - The engine converts 97.7 % of furfural to FFT within ~15 s and destroys ~99 % of it.
  - The cysteine-over-H2S FFT contrast is predicted 0.234 vs 0.046 measured (+0.71 dex).
  - The H2S dose order passes (−0.03 dex).
  - Source: CF-065, Munich cluster.
- **Sig:**
  - `k_fur_decay` and `k_fur_fft` clipped on the upper side in the B60 profile;
  - `k_fur_fft_hs` is unpinned;
  - right level, wrong turnover.
- **Step 4:**
  - Integrity R: Munch T5 and Hofmann T10 read on the page.
  - Conditions: speciation and headspace remove ~1.1 dex of the Zhai gap.
  - Mechanism: the contrast constrains furfural lifetime × `k_cys_h2s` × FFT survival jointly.
- **Cand:**
  - (a) H2S vessel partition is missing (Vg/Vl 3: 15–24 % dissolved; `k_fur_fft_hs` is biased low by 0.6–0.8 dex);
  - (b) the furfural lifetime;
  - (c) FFT survival in the hold.
- **Shared:** `r_fur_fft_hs`, `r_fur_fft`, `k_fur_decay`, `k_fft_decay`; G9 (headspace), G10.
- **Checks:** CK-34, CK-35, CK-36.
- **W5.4:** Th × cys, Th × ck.
- **Decisive:** the CF-065 deciders: (ii) residual furfural and H2S at 20 min; (iii) direct FFT stability.

**MF-14 — Absolute thiol levels for a new laboratory (H0)**

- **Obs:** MFT and FFT levels across laboratories.
- **Size:**
  - τ_ABS 0.964 dex, U ≥ 1.60 (`kinetic_core_thiol_lolo_nested_tau.md`).
  - Yiltirak hold-out MFT +2.44, FFT +1.98 dex.
  - Yiltirak vs Schieberle x120 (CF-068).
  - YIL-01: the engine rises with T where the measurement falls (−0.60 dex).
- **Sig:** leave-one-group-out tension; one sign on Yiltirak.
- **Step 5:**
  - Applied as a refusal by name (EX-059, D2 branch (i)).
  - Integrity R: both laboratories' units verified by two readers.
  - Conditions: engine-free arithmetic explains 0.35–0.95 of the 2.12 dex gap.
  - The O2-residual regression has R² 0.56, n 4.
- **Cand:**
  - (a) headspace O2 (17 mL air over 3 mL);
  - (b) post-heating thiol loss in Yiltirak's work-up;
  - (c) DFA reads high (P ≥ 0.5 dex: open).
- **Shared:** G9 (a vessel / headspace operator), G3.
- **Checks:** CK-37, CK-38.
- **Blind:** the WUR blind-test samples are BLIND; never read them for this.
- **W5.4:** Th × every lever (absolute levels are refused; only ratios ship).
- **Decisive:** DM:thiol_absolute_second_lab.

**MF-15 — MFT barely responds to the cysteine / GSH dose in thiamine pots**

- **Obs:** Grosch & Zeiler-Hilgart 1992 Tables II / III (thiamine + Cys, SIDA); Zhou 2023c GSH; Bolton 1994.
- **Size:**
  - Grosch Cys no1 vs no4: measured −0.85 (x7), predicted −0.00.
  - Three contrasts have the wrong sign.
  - Z23C-GSH-01 −1.12 dex, wrong sign.
  - The Bolton 1994 thiamine pot is +1.30 dex (232.9 vs 11.7 ppb; case M4).
  - Source: `kinetic_core_contrast_benchmark_score_v1.md`; CF-089.
- **Sig:** identical predictions on the moved axis (no term); wrong sign.
- **Step 2:**
  - Integrity K (Grosch is a keyed held-out "clean" table).
  - Conditions: CF-089 / CF-090 (dose, time, medium, pH).
  - Never checked: R for the engine's pot definitions (thiamine and Cys charges).
- **Cand:**
  - (a) the thiamine route is first order in thiamine with no cysteine term (CF-089 remedy);
  - (b) instant Cys supply (MF-02);
  - (c) GSH has only `r_gsh_h2s` (F27-c).
- **Shared:** "THE THIAMINE ROUTE" (sulfur.py), `r_gsh_h2s`; G4.
- **Checks:** CK-10, CK-42.
- **W5.4:** Th × thi, Th × cys.
- **Decisive:** a cysteine × thiamine series with thiamine measured (CF-089). No DM entry.

**MF-16 — Disulfide share a decade short; `k_dimer_mft` on its instant bound**

- **Obs:** bis(2-methyl-3-furyl) disulfide / MFT (Zhou 2023, Zhang 2024 Cys arm).
- **Size:**
  - x11 / x8.8 short at B9 (≈ −1.0 dex);
  - at B60, −2.71σ (Zhou 2023) and −1.99σ (Zhang 2024).
  - Source: case M1.
- **Sig:** `k_dimer_mft` is 0.49999, ACTIVE on its upper bound 0.5 (B60 profile); one sign in two laboratories.
- **Step 3:**
  - Measurement: the M1 review.
  - Never checked: R (the species definitions and the basis of the share).
- **Cand:**
  - (a) the oxidant pool limits (the deoxypentosones as oxidant; C3);
  - (b) the formation rate is not the limit.
- **Shared:** `k_dimer_*`; G3.
- **Checks:** CK-09, CK-28.
- **W5.4:** Th × cys, Th × ck.
- **Decisive:** DM:thiol_supply_structure (C3 side); atlas E12.

**MF-17 — 2-Acetylthiazole far too low, no route, wrong pH sign**

- **Obs:** ACTZ, Zhang 2025b and Xu 2007.
- **Size:**
  - −3.25 dex (114 vs 0.065; `core_panel_scores.md`).
  - The R3En offset on `k_cys_actz` is −2.64.
  - Xu 2007 pH 5.6 / 7.0: measured −0.26, predicted +2.31 (wrong sign).
- **Sig:** one sign; wrong sign with pH.
- **Step 1:**
  - Integrity is partial: CF-009 (a possible double count by column sum, ≤ 0.3 dex).
  - Hold-out status is inconsistent: the Zhang 2025b blank volatiles are a frozen hold-out (Amendment 45), yet the
    ACTZ row is on the panel.
  - Never checked: R on the hold-out status; conditions.
- **Cand:**
  - (a) the thiazoline → thiazole route (R37–R39, listed, not modelled), which needs NH3 (G1 amine fate);
  - (b) thiazoles rise with pH (P15).
- **Shared:** `r_cys_actz`, the amine fate (G1).
- **Checks:** CK-07, CK-22.
- **W5.4:** Th × pH, Th × cys.
- **Decisive:** none written.

**MF-18 — Lipid is inert on the thiols**

- **Obs:** MFT and FFT with lipid (Whitfield 1988, Elmore 2002, Xu 2011, Yiltirak 2026, US 8,920,862).
- **Size:** the engine predicts 0, against measured falls of x3–143 (−0.5 to −2.2 dex). A partition term explains at
  most 0.1–0.4 dex (CF-031; case M5).
- **Sig:** identical predictions on the moved axis.
- **Step 4:**
  - Measurement: species and aliquot (resolve2_sulfur; the pre-registered speciation rule FAILS).
  - Mechanism: electrophile capture is proposed at P 0.55, not applied.
- **Cand:**
  - (a) lipid-electrophile capture of the MFT pool;
  - (b) oil-phase partition (bounded, ≤ 0.4 dex);
  - (c) aliquot sampling in Yiltirak.
- **Shared:** lanes summed directly (lipid.py); the thiol-sink law (G3).
- **Checks:** CK-29.
- **Reserved:**
  - the E7 TEST contrasts F30#11 (Kim 1996b), F30#15 (Zamora 2015) and F30#16 (Hidalgo 1993) are UNSEEN;
  - Xu 2011 is E12 TEST (item 7 only);
  - `heldout_zhang2021_table3.yml` is keyed held-out.
- **W5.4:** Th × fat (the whole fat column for thiols).
- **Decisive:** atlas E7.

### Trunk lane

**MF-19 — Pentose + amino-acid sugar steps carry no pH**

- **Obs:** NF and furfural in the Meynier 1995 Gly / Met + ribose ladder; xylose + Gly flags in R3En.
- **Size:**
  - Meynier NF measured −0.02 to +1.21 dex, predicted 0.00 on every pair; FUR up to +2.70 against 0.00.
  - R3En: desclaux2006 +1.41σ (8 rows) and davidek2008 +1.36σ (19 rows).
  - Source: `kinetic_core_contrast_benchmark_score_v1.md`.
- **Sig:** identical predictions on the moved axis; correlates with pH.
- **Step 4:**
  - Mechanism: case M1b (rivals R0–R5); E2 is WEAKENED (P 0.55).
  - Conditions: CF-118 (the Martins slopes are on the label pH; 0.69 → ~0.55 on the mean pH).
- **Cand:**
  - (a) `TRUNK_PH_TERMS` with a measured slope on the mean pH;
  - (b) the split moves with pH vs 3-DX removal (R2 vs R5).
- **Shared:** G2; G1.
- **Checks:** CK-19 (owned by M1b).
- **Reserved:** Meynier 1995 (`heldout_meynier1995_tab1_2_4.yml`) is read only in M1b item 7.
- **W5.4:** Th × pH, Pz × pH, St × pH.
- **Decisive:** DM:pentose_split_vs_removal.

**MF-20 — Deoxyosone partition: DPO low, TDP high, the split's pH sign wrong**

- **Obs:** 3-DX / 1-DX in Jiangnan pots; Zhai 2023 hold-outs.
- **Size:**
  - R3En: DPO −0.37σ over 112 rows (one sign); TDP +0.42σ over 127 rows.
  - Zhai 2023 JAFC DPO −1.52 dex; Zhai 2023 FoodChem −1.16 dex.
  - The engine 3-DX / 1-DX rises with pH where the Jiangnan data fall (EX-033).
- **Sig:** one sign across many rows; correlates with pH.
- **Step 5:**
  - The split law is identified (B7S), but its pH slope stays ~0 (EX-036, EX-046, EX-053).
  - CF-071: signature (iii) is not thiol-diagnostic.
  - Decisive measurement only; no seen-data discriminator (DM entry: "none is on disk").
- **Cand:** R2 (the split moves) vs R5 (base-favoured 3-DX removal).
- **Shared:** G2.
- **Checks:** CK-20 (none runnable on seen data).
- **W5.4:** Th × pH.
- **Decisive:** DM:pentose_split_vs_removal.

**MF-21 — Xylose-alone dicarbonyls at −16σ; the laboratory's dicarbonyl scale spans 10x**

- **Obs:** zhou2025b xylose-alone dicarbonyls (R3En); CF-060 / CF-082 AX-ARP pots.
- **Size:**
  - −16.12σ over 4 rows (R3En). The size suggests a zero or floor prediction.
  - CF-060: the five AX-ARP pots span ~10x.
- **Sig:** one sign; implausible size.
- **Step 1:**
  - Integrity R for CF-060 (nitrogen closure condemns the two HIGH prints, P 0.6).
  - Never checked: R for the −16σ rows themselves (the engine output and the keyed rows).
- **Cand:**
  - (a) the engine predicts ~0 (a pot-definition or route gap);
  - (b) the scale issue of CF-060.
- **Shared:** the trunk xylose-alone pathway.
- **Checks:** CK-02.
- **W5.4:** Th × sug (indirect).
- **Decisive:** the OPD dilution factor of each paper (CF-060).

**MF-22 — Acetic-acid loss ≥ 40x too fast**

- **Obs:** acetic acid (Martins 2005 fit; Davidek 2008 plateau).
- **Size:**
  - `r_aa_frag` 0.0164 /min at 120 C, against ≤ 4e-4 measured (≥ 1.6 dex).
  - R3En AA +1.62σ.
  - Its Ea sits on the lower bound of 20 kJ/mol.
  - Source: CF-069.
- **Sig:** constant on a bound; one sign.
- **Step 4:** cluster 3 arm B (pH-free trunk; with `k_aa_frag` off and inline pH it fits RMS 0.16 and predicts
  Martins 2005b held-out). Explained for the AA rows, P 0.6.
- **Cand:**
  - (a) the missing pH state (G1);
  - (b) a late acid sink or a formic / acetic split.
- **Shared:** G1.
- **Checks:** CK-16 (owned by M1b, D6 Davidek).
- **W5.4:** none directly. Acid budget → pH × all.
- **Decisive:** arm E fit (CF-069).

**MF-23 — Glucose alone in an autoclave: osones and glyoxal low (one pot)**

- **Obs:** Steinhagen 2021, 121 C.
- **Size:**
  - glucosone −1.80, glyoxal −1.54, HMF −0.97, 3,4-DDG −0.82 dex;
  - 3-DG +0.10; MGO 0.00.
  - TEMP-05 and SUG-12 have the wrong sign.
  - Source: `core_panel_scores.md`.
- **Sig:** one sign on the oxidative products only (a pattern, one pot).
- **Step 1:** never re-read (R).
- **Cand:**
  - (a) no oxidative osone route (`ATMOSPHERE_FACTORS` = {air: 1.0}; F08-a);
  - (b) the amine-free lane needs phosphate (F14-a).
- **Shared:** G9 (O2), G8 (osone supply).
- **Checks:** CK-03.
- **W5.4:** St × sug, Pz × sug (donor supply).
- **Decisive:** none written.

**MF-24 — MGO self-loss: first order with Ea on its floor vs second order at 95–121 kJ/mol**

- **Obs:** MGO loss (trunk `k_mgo_mel` vs Kanzler 2026).
- **Size:** `k_mgo_mel` gives 0.0158 /min at 80 C with Ea 20.04 kJ/mol. Kanzler's apparent Ea is 95–121 kJ/mol
  (MGO + amino acid), above even B1's CI. Amino acids speed MGO loss 2.6–9.7x (CF-122). The fold in k at a given
  temperature is not computed in any record.
- **Sig:** constant on a bound (FITTED_EA_BOUNDS floor; B1 CI −42 to 82).
- **Step 4:** CF-122 (concentration: Kanzler at 100 mM, the trunk ~1000x more dilute).
- **Cand:**
  - (a) `k_mgo_mel` lumps MGO + amine chemistry;
  - (b) the order differs with concentration.
- **Shared:** `r_mgo_mel`; G6 (adducts).
- **Checks:** CK-46, CK-47.
- **W5.4:** St × Tt, Pz × Tt.
- **Decisive:** DM:methylglyoxal_self_loss.

**MF-25 — Amino-acid consumption by fed MGO: the engine is ≥ 1e3x low**

- **Obs:** glycine loss in Piloty & Baltes 1979 (MGO 312 mM + 13 amino acids, 100 C, 3 h).
- **Size:** 82–86 % measured vs 1.9e-5 to 1.6e-3 in the engine: ≥ 3 dex (a lower bound; CF-095). Sedehi 2013's gap is
  ~2e4.
- **Sig:** identical across pH 4.5–9.5 (a pH-flat consumption the engine lacks); two laboratories.
- **Step 4:** resolve2_pyrbuf Q2-b (units red-teamed).
- **Cand:**
  - (a) the missing dicarbonyl-amine adducts (imidazoles, hydroimidazolones);
  - (b) Piloty's mixture and buffer inflate consumption (a factor, not 1e3).
- **Shared:** G6.
- **Checks:** CK-46.
- **W5.4:** St × cys, Pz × cys (the amine budget), St × sug.
- **Decisive:** none written. The candidate A3 step uses Piloty's shape check.

**MF-26 — No buffer catalysis (E5)**

- **Obs:**
  - buffered vs unbuffered MFT (Mottram 2002: x2.2–2.6, FFT unchanged);
  - DFG decay (Davidek 2002: x2.8–3.6);
  - ARP decay (Deng 2021, CF-120).
- **Size:** 0.34–0.56 dex on the moved axis; the engine gives 0.
- **Sig:**
  - identical predictions on the moved axis;
  - step-specific (FFT and acetic acid unchanged);
  - 8+ laboratories.
- **Step 4:** CF-093 (non-saturating to 0.5 M in dilute pots, P 0.8); CF-094 (carboxylate selectivity, P 0.5);
  CF-120.
- **Cand:**
  - (a) a per-step phosphate factor on resolved steps only;
  - (b) general catalysis on every step (rejected for saturation; selectivity open).
- **Shared:** G5.
- **Checks:** CK-54 (owned by M5), CK-55 (done).
- **Reserved:**
  - the E5 TEST set (Potman, Kallen 1971b, Thornalley 1999) is already SPENT (EX-058);
  - Larralde 1995 is keyed held-out.
- **W5.4:** every odorant × ck (buffer and matrix of the food).
- **Decisive:** CF-011 buffer arm C with a citrate arm; atlas E5.

**MF-27 — Protein glycation: barrier too low; no lysine plateau**

- **Obs:** FDNB available lysine in soy (Wolf 1977a and its companion).
- **Size:**
  - E 137.7 vs the engine's 96.8 kJ/mol (~5 SE).
  - The engine over-predicts by +0.24 to +0.42 dex at 100 C and +0.5 to +0.7 dex at 80 C.
  - On pass 2 it loses another 36 %, against ~0 measured.
  - Source: CF-100; CF-096.
- **Sig:** far from its measured analogue; correlates with temperature; wrong shape (plateau).
- **Step 4:** Wolf 1977a findings, red-teamed.
- **Cand:**
  - (a) plateau-capped end points inflate Wolf's E;
  - (b) an inaccessible lysine fraction;
  - (c) FDNB counts bound Amadori lysine.
- **Shared:** `k_glyc`, `r_flp_decay` (B20).
- **Checks:** CK-56.
- **W5.4:** none directly (protein amine supply → Pz, St × ck).
- **Decisive:** a soy or pea glycation time course at ≥ 3 temperatures (CF-100).

### Strecker lane

**MF-29 — Argon PA + PAA: every M3 rival ~5x high (no candidate on record)**

- **Obs:** Hofmann & Schieberle 2000b Table 2, expt I. Fed ARP-Phe 100 mM, 100 C, 120 min, closed vial, argon.
- **Size:**
  - Measured total 0.09 mM (PA 0.06 + PAA 0.03; read error 0.036 / 0.072 dex).
  - Median prediction +0.63 to +0.73 dex in all four rivals: R1 0.468, R2 0.465, R3 0.479, R4 0.387 mM.
  - The q10 reaches the datum only in R2 (+0.07).
  - Source: `kinetic_core_m3_phase1_checks.md` on `m3-phase1-build`, "Open diagnostic: Table 2's argon PA + PAA";
    prereg v2.1 DRAFT sec. 22 A6.
- **Sig:**
  - one sign across every rival;
  - a free prediction (σ_ox is solved to the AIR EXCESS, so the argon total is the anoxic part alone: the e1 / e3 / eF
    exits feeding route A).
- **Step 1:**
  - Integrity K only (two-reader, adjudicated 2026-10-31). **R never done.**
  - Conditions are recorded but never related to this misfit:
    - argon deoxygenation and headspace are not printed;
    - the air arm reads "flushing with oxygen" (≈ 4.8x O2 if pure O2);
    - phosphate is printed 0.5 mmol/L in the methods vs 0.5 mol/L in the footnote (CF-020);
    - ARP purity is ~90 %.
  - The power check makes argon PA SYNTHETIC from the truth, so a common 0.7 dex miss is invisible to power.
  - Note: U3a, the primary claim block, divides by this argon PA.
- **Cand:** none named by any record for this misfit. Adjacent open items in the M3 prereg, not named as its causes:
  - (a) the k_d band vs the slow fed decay (CF-121: 0.006–0.010 vs ≥ 0.0126 /min);
  - (b) the s_DG1 upper end (Granvogl x10^0.82, flagged "generous to R2");
  - (c) a donor sink κ_D[D][Phe] (D-b) or the adducts of MF-25;
  - (d) the air arm being O2-flushed changes σ_ox, and so k_anox = k_d (1 − σ_ox);
  - (e) ARP purity (≤ 0.05 dex).
- **Shared:** k_d, k_anox, the anoxic shares a1 / a3 / aF, the route-A shares (s_DG1, s_FOX, s_3DG, calibrated on
  hofmann2000 Table 4 argon); G6, G9, G11.
- **Checks:** CK-01, CK-43, CK-44, CK-50.
- **Spent:** Table 2 is SEEN in full (B29, EX-049).
- **Blind:** hofmann2000b Figs. 2–3 (trapped osones and Phe, argon and air + Cu, 85 C) are BLIND to M3. No check here
  may read them.
- **W5.4:** St × ck (atmosphere and headspace), St × sug, St × Tt.
- **Decisive:** DM:strecker_routes (a) (fed ARP-Phe, argon vs air, headspace stated). The entry does not mention this
  misfit.

**MF-30 — The B18 Strecker constant is ~2.3 dex below Hofmann's fed MGO + Phe; the glucose-pot PA has a lag shape**

- **Obs:** PA in the glucose + Phe pot (hofmann2000 Fig. 2) under the R0-pinned (B18) constant; fed MGO + Phe
  (EX-021).
- **Size:**
  - −3.04 / −2.27 / −2.01 / −1.88 / −1.83 / −1.70 dex at 10–500 min (`kinetic_core_m3_phase1_screen.md` on
    `m3-phase1-build`, by arithmetic);
  - EX-021 2.3 dex.
- **Sig:**
  - constant far from its measured analogue;
  - wrong shape (the deficit shrinks with time);
  - E4's step pH law double counts (k5 / k7 0.069).
- **Step 4:**
  - EX-021.
  - Case M3 lists compensators: B18's constant was fitted on pyrazine-pot rows.
- **Cand:**
  - (a) the step constant (fitted on B18 pyrazine rows, absorbing the supply deficit of MF-34);
  - (b) the donor supply (G8).
- **Shared:** `r_go_ak`, `r_mgo_ak`, PYRAZINE_PH_SLOPES; G8.
- **Checks:** CK-45.
- **W5.4:** St × every lever; Pz × pH.
- **Decisive:** DM:strecker_routes.

**MF-31 — Opposite signs across laboratories (Pan low, Deng high); peak time wrong**

- **Obs:** Strecker aldehydes and methional (B22, B22b).
- **Size:**
  - Pan 2025 −3.6 to −5.6 dex; Deng 2022 +1.5 to +2.9 dex (methional +1.93 at 120 min).
  - The model peaks at 30 min against 120 measured.
  - Source: case M3.
- **Sig:** leave-one-group-out tension (opposite signs); wrong shape.
- **Step 1:** never re-read (R). Pan's 3.6–5.6 dex is the size of a unit or basis slip.
- **Cand:**
  - (a) a unit / basis issue (Pan);
  - (b) amine identity and methionine handling (`METHIONINE_SHIPPED = False`);
  - (c) the donor supply.
- **Shared:** `r_mgo_ak`, the methionine Strecker step.
- **Checks:** CK-08 (Pan half only), CK-48.
- **Reserved:** Deng 2022 is M3 phase-2 reserved (U2-Deng, "reported, engine-seen"). The Deng half is NOT RUNNABLE
  until phase 2 freezes.
- **W5.4:** St × cys (methional), St × Tt.
- **Decisive:** DM:strecker_routes.

**MF-32 — ARP-Phe decay: the fed and in-pot rates disagree; the fed early / late shape drifts**

- **Obs:**
  - ARP-Phe in the glucose + Phe pot (Table 3);
  - fed ARP-Phe PA shape (hofmann2000b Fig. 1; EX-051);
  - phosphate on ARP decay (CF-120).
- **Size:**
  - In-pot k_d ≥ 0.0126–0.0144 /min, against fed ~0.006–0.010 (0.1–0.35 dex).
  - EX-051 G7: +0.5 dex at 10 min and −0.4 to −0.5 dex at 500 min.
  - The PA / PAA drift is 3.33 → 1.11.
- **Sig:** leave-one-group-out tension inside one laboratory; right level, wrong shape.
- **Step 3:**
  - Integrity R: CF-121 (two fed designs, not three; the 0.5 mmol/L print).
  - Measurement: at 10 % analyser CV the tension weakens.
- **Cand:**
  - (a) an order or vessel effect (O2-limited reflux in Fig. 1);
  - (b) the acid load in the trapped pot (0.07–0.7 dex);
  - (c) noise.
- **Shared:** k_d, ψ; G1 (acid load), G11.
- **Checks:** CK-49.
- **Blind:** Figs. 2–3 would decide; NOT RUNNABLE.
- **W5.4:** St × Tt, St × ck.
- **Decisive:** DM:strecker_routes (a).

**MF-33 — O2 dependence of Strecker output: the engine air / argon ratio is 2.19 against 9.2; the R2 layer needs more O2 than the vial holds**

- **Obs:** Hofmann 2000b Table 2 air / argon (B29); the M3 R2 calibration.
- **Size:**
  - The engine's air / argon is 2.19 against 9.2 (−0.62 dex); the branching order is wrong (0.79x vs 2.63x).
  - R2 is feasible in 66 / 600 draws; `c_Cu5` is at its bound 25 in 362 draws.
  - R2 needs 250–700 µmol O2 in the air arm, against R1's 4–8.
  - Source: case M3; `kinetic_core_m3_phase1_screen.md`.
- **Sig:** constant on a bound; wrong ratio.
- **Step 2:** the headspace is unprinted (M3 open question 11).
- **Cand:**
  - (a) the air arm is O2-flushed (≈ 4.8x O2);
  - (b) R2's glucosone route is wrong;
  - (c) a third O2 route (Wietstock / Voigt, EX-067).
- **Shared:** G9.
- **Checks:** CK-50, CK-01.
- **W5.4:** St × ck.
- **Decisive:** DM:strecker_routes (a).

### Pyrazine lane

**MF-34 — Total pyrazines ~2.9 dex low: a supply deficit**

- **Obs:** pyrazine totals (B18 T4, B21 Leahy).
- **Size:** −2.9 dex (T4); −2.8 dex (B21). The Xia trunk MGO is 7.8 mM against ≤ 0.5 mM (+1.19 dex, opposite sign).
- **Sig:** one sign; glycine stands in for lysine (amine identity).
- **Step 4:** case M6 ("a dicarbonyl-supply and amine-identity deficit, not a sink").
- **Cand:**
  - (a) dicarbonyl supply;
  - (b) amine identity (Lys vs Gly);
  - (c) the aminoketone ring step (R1 vs R3).
- **Shared:** G8, G6 (adducts would deepen the deficit: a sign check).
- **Checks:** CK-51 (owned by M6).
- **Reserved:**
  - the E13 TEST studies (Shibamoto 1977, Low 2007, Van Lancker 2012, Amrani-Hemaimi 1995) are owed and UNSCORED;
  - the E14 TEST set;
  - Lee 2024 is keyed held-out.
- **W5.4:** Pz × every lever.
- **Decisive:** DM:pyrazine_ring_aminoketone.

**MF-35 — Pyrazine : 2,5-DMP 23 : 1 measured vs 6e-8 : 1; the mixing ratio R**

- **Obs:** product ratios (Leahy, B18 T3; Zhou 2023c glycine column, CF-087; UC Davis, CF-098).
- **Size:**
  - −8.6 dex on the ratio.
  - `k_cond` is 1e4 with Ea 0, against an analogue ~1e7x slower with Ea 44 (E14).
  - Leahy rows are +0.51 dex under CF-077.
- **Sig:** far from its measured analogue; a forbidden low R with glycine.
- **Step 4:** EX-052, EX-057; CF-087 (low R with glycine is a condition effect).
- **Cand:**
  - (a) `k_cond` at the analogue's measured constant;
  - (b) segregated formation (h < 0.05);
  - (c) low-level data.
- **Shared:** `r_ak_mpz`, `r_akg_pz`, `r_akm_dmp`, `k_cond`.
- **Checks:** CK-52.
- **Reserved:** as MF-34. UC Davis is an E13 TEST laboratory: read only as already quoted.
- **W5.4:** Pz × sug, Pz × cys (amine).
- **Decisive:** atlas E13–E14.

**MF-36 — The pyrazine pH law double counts on the Strecker step**

- **Obs:** methylpyrazine pH 7 / 5 (Leahy).
- **Size:** shipped 170 vs measured 14 (+1.08 dex).
- **Sig:** correlates with pH.
- **Step 4:** atlas sec. 5. With the step slope at 0.07, methylpyrazine matches (13.8) but pyrazine goes flat (0.71).
- **Cand:**
  - (a) the slope belongs to the supply only;
  - (b) the step slope is real, below pH 7 (E4 design).
- **Shared:** PYRAZINE_PH_SLOPES; G2.
- **Checks:** CK-53 (done; the remaining pyrazine-flat side is in M6).
- **W5.4:** Pz × pH, St × pH.
- **Decisive:** atlas E4.

### Acrylamide lane (incl. furanics)

**MF-28 — Furanics are pH-insensitive and low (glucose + Ala, 130 C)**

- **Obs:** DMHF and HMF, Schibilsky 2019, pH 5 vs 8.
- **Size:**
  - DMHF −1.72 dex (pH 5) and −2.43 dex (pH 8): identical predictions (21.841) against a measured rise of x5.1.
  - HMF: identical at both pH values, against a measured x1.76.
  - Source: `core_panel_scores.md` (lane column "acrylamide").
- **Sig:** identical predictions on the moved axis.
- **Step 1:** never re-read (R). Identical output at two pH values may be the pot not passing pH to the lane (check
  first).
- **Cand:**
  - (a) no pH term in the lane;
  - (b) 2,3-enolisation rises with pH (E2);
  - (c) HDMF phosphate x4.3 (E5).
- **Shared:** G2, G5.
- **Checks:** CK-06, CK-21.
- **W5.4:** AA × pH (furanics share the acrylamide carbonyl supply).
- **Decisive:** none written.

**MF-37 — Glucose + Asn, 180 C: too high early, wrong time shape**

- **Obs:** acrylamide, Chang 2021 (10 / 30 min, water and buffer); Ye 2024.
- **Size:**
  - Chang 10 min +2.34 dex; 30 min +0.40; 30 min water +0.65.
  - Measured 10 → 30 min is x52; the engine falls.
  - Ye 2024 +1.70.
  - TEMP-01 has the wrong direction.
  - The contrasts: 0 / 5 signs right.
- **Sig:** right level late, wrong shape; correlates with temperature.
- **Step 1:** never re-read (R) for the heat-up, the units or the basis.
- **Cand:**
  - (a) heat-up / come-up time (conditions);
  - (b) the phosphate pH law applied to every buffer (`acrylamide_conditions.py`);
  - (c) irreversible amine adducts.
- **Shared:** G5 (buffer law).
- **Checks:** CK-04, CK-57.
- **W5.4:** AA × Tt, AA × ck.
- **Decisive:** none written.

**MF-38 — Fructose + Asn reads low**

- **Obs:** acrylamide and HMF, Lin 2022.
- **Size:** acrylamide −0.96, HMF −0.80 dex.
- **Sig:** one sign on the sugar axis. Asn reacts only with Glc in the lane, against fructose > glucose x1.1–2.9 in 3
  laboratories (F29-a).
- **Step 1:** never re-read (R).
- **Cand:** (a) a fructose + Asn step at the measured Fru / Glc ratio.
- **Shared:** `a_asn_glc_sb`.
- **Checks:** CK-04, CK-58.
- **W5.4:** AA × sug.
- **Decisive:** atlas E9.

**MF-39 — Extruded soy protein isolate at 130 C: −3.98 dex**

- **Obs:** acrylamide, a trust_loop "predictive" row (ACSRef3).
- **Size:** −3.98 dex (150 measured vs 0.016).
- **Sig:** single row, very large: an integrity signature.
- **Step 1:** none recorded.
- **Cand:**
  - (a) the precursor charge (the SPI's Asn and sugar);
  - (b) a units slip;
  - (c) matrix / a_w.
- **Shared:** the lane's precursor mapping.
- **Checks:** CK-04, CK-59.
- **W5.4:** AA × ck (extrusion).
- **Decisive:** none written.

**MF-40 — The acrylamide lane's emergent Ea at 90–120 C: 296 vs 103 kJ/mol**

- **Obs:** low-temperature yields (Gokmen 2012; Knol 2005 EaF 94.4, predicted 189.5).
- **Size:** 60-min yields x2000 low at 90 C (−3.3 dex) and x15 low at 120 C (−1.2 dex).
- **Sig:** correlates with temperature; outside the fit window (120–200 C).
- **Step 5:** the guard is applied (`LOW_T_REFUSAL_C`); CF-039 estimand work.
- **Cand:**
  - (a) a second, lower-barrier route below ~120 C;
  - (b) Gokmen's estimator inflates low-T yields (partly).
- **Shared:** `k_int1_acr`, `k_asn_glc`.
- **Checks:** CK-60 (none discriminating on seen data).
- **W5.4:** AA × Tt (low-temperature cooking).
- **Decisive:** acrylamide + Asn time courses at ≥ 3 temperatures in 80–120 C (CF-106).

**MF-41 — Water activity: the wet arm is 2.6x slow; elimination falls dry-side where the data rise x6–14**

- **Obs:** k_F and k_E, Bassama 2011 (plantain, a_w 0.97 / 0.90 / 0.43).
- **Size:** formation wet arm −0.41 dex; elimination residual x3.1 / x6.7 (0.5–0.8 dex).
- **Sig:** correlates with a_w.
- **Step 4:** CF-104 (wet-basis mass action supported for the dry arms; the flat multiplier is upheld at P 0.75).
- **Cand:**
  - (a) the moisture range is outside calibration;
  - (b) a lumped k_F;
  - (c) elimination not first order.
- **Shared:** `AW_MULTIPLIER`, the elimination a_w shape.
- **Checks:** CK-61 (done; no seen-data discriminator left).
- **W5.4:** AA × ck (crust moisture), AA × Tt.
- **Decisive:** none written.

### Lipid lane

**MF-42 — Lipid aldehydes and 2-pentylfuran read low in matrices**

- **Obs:**
  - li_2026 SPI–WG HME: 2-pentylfuran −2.56, hexanal −0.94, nonanal −0.49 dex;
  - bi_2020 roasted pea: hexanal −0.47;
  - Trikusuma: −0.34 / −0.40 / −0.19 / +0.21.
  - Source: `core_panel_scores.md`.
- **Size:** lane median −0.45 dex, sign-consistent in 88 % of rows (`lane_offset_diagnostic.md`, stale, B46).
- **Sig:** one sign. Only 3 pots and 2 temperature levels (`data_wishlist.md` sec. 7).
- **Step 1:** never re-read (R).
- **Cand:**
  - (a) the hydroperoxide / precursor charge of the matrix;
  - (b) matrix binding (K_g);
  - (c) extraction conditions (cf. CF-022, a headspace-condition error found before in this lane).
- **Shared:** the matrix layer (`HeadspaceModel`); the lipid precursor charge.
- **Checks:** CK-05, CK-62.
- **W5.4:** La × fat, La × ck.
- **Decisive:** none written.

(MF-28 is listed under acrylamide because the panel scores DMHF on that lane.)

---

## 3. Every leading candidate's cheapest check on seen data, in batches

Each check states its own unit, noise band and in-sample rows; the pass rule is the one in the header. *Work* is
wall-clock hours of work and *cpu* is unattended single-core compute.
- *Spends* lists the datasets whose values the check reads. Datasets already SEEN are re-spent at no new loss, but
  must still be declared.
- **owned by X**: the check reads a case's decisive data and runs only inside that case.
- **done**: a record already holds the check's result. It is cited, not rerun.

The batches are ordered by shared structure, cheapest discriminating first. Inside a batch, run the checks in the
order listed.

### Batch 0 — integrity (step 1) for the misfits never re-read

No model run; PDFs are in the main checkout's `data/articles/` (laptop). These checks re-read sources whose values are
already seen. Where a source is a held-out bundle, read its methods and the scored cells only. Total ≈ 9.25 h.

| id | MF | what to read / compute | work | spends | pass rule specifics |
|---|---|---|---|---|---|
| CK-01 | 29, 33 | Hofmann 2000b Table 2 and Methods. Check: <ul><li>argon total vs PA / PAA split;</li><li>µmol per mmol of ARP or of Phe;</li><li>vial volume and fill;</li><li>deoxygenation method;</li><li>"flushing with oxygen" vs air;</li><li>0.5 mol/L vs mmol/L phosphate;</li><li>ARP purity.</li></ul>Compare with `strecker_hofmann2000b.yml` and the dossier | 0.75 | Table 2, Methods (SEEN) | Integrity fails if any keyed cell, basis or condition differs from the page. If the argon cell's basis differs, MF-29 is explained if the corrected datum lies within the rivals' q10–q90 |
| CK-02 | 21 | R3En's zhou2025b xylose-alone rows vs the PDF; the engine's prediction for that pot (is it ~0 or a floor?) | 0.75 | zhou2025b (SEEN) | Explained (integrity / pot definition) if the −16σ comes from a zero prediction or a mis-keyed unit |
| CK-03 | 23 | Steinhagen 2021: autoclave profile, headspace / O2, pH, units of glucosone and glyoxal | 0.75 | Steinhagen (SEEN, maillard_path_holdout) | A condition the pot omits (O2, heat-up) is a step-2 lead; a unit slip of ≥ 1 dex explains |
| CK-04 | 37, 38, 39 | Chang 2021 (heat-up, buffer, 10 / 30 min), Ye 2024 (µmol/mol basis), Lin 2022, the extruded-SPI row (units, precursor charge) | 2.0 | the four rows' sources (SEEN) | As CK-03, for each row |
| CK-05 | 42 | li_2026, bi_2020, Trikusuma: extraction, headspace, matrix, units | 1.5 | the three sources (SEEN) | As CK-03 |
| CK-06 | 28 | Schibilsky 2019 bundle: is the pH passed to the lane at all (identical outputs at pH 5 and 8)? Then the units and basis | 0.75 | Schibilsky (SEEN) | Explained (pot-definition bug) if the pH is not passed; otherwise the misfit stands as structural |
| CK-07 | 17 | `FIT_HOLDOUT_DECLARATION.md` Amendment 45 vs the panel's Zhang 2025b ACTZ row: is the row a frozen hold-out that is being scored? CF-009 column sum | 0.5 | none new (declaration only) | Records the status; flags a hold-out breach if one exists |
| CK-08 | 31, 35 | **Pan 2025 B22 rows** (units, basis); Leahy T3 rows (the B18 fit rows) | 1.0 | Pan 2025, Leahy (SEEN) | A ≥ 3 dex unit slip explains Pan's sign. The Deng 2022 half is NOT RUNNABLE (M3 phase 2 reserved) |
| CK-09 | 16 | Zhou 2023 / Zhang 2024 disulfide rows: the species (bis(2-methyl-3-furyl) disulfide vs total), the share's basis | 0.75 | both (SEEN) | A basis mismatch of ≥ 0.5 dex explains part |
| CK-10 | 15 | Grosch 1992 Tables II / III vs the engine's pot definitions (thiamine and Cys charges, pyrophosphate, reflux) | 0.5 | Grosch T2 / T3 (SEEN, scored) | Explained if the pots carry the wrong charges |

### Batch A — the in-situ pH state (`ph_state`; G1)

One forward machinery serves MF-06, -08, -03, -22 and -32. Total ≈ 5.25 h runnable work plus 1.5 h owned by M1b,
< 0.5 h cpu.

| id | MF / cand | what to compute | work | spends | unit · noise band · in-sample rows |
|---|---|---|---|---|---|
| CK-12 | 06 (a); side readings 03, 08, 17 | **Owned by M1b (coordinator, 2026-11-06): it spends Hou and Meynier, CF-084's controls and so M1b(i)'s decisive data; run it only inside M1b.** Forward cooled-pH, no refit. <ul><li>The amine is returned as NH3 on `r_cys_h2s` / `r_cys_thermal` / `ch_cys_ox`, at Sohn 1995's measured release.</li><li>The ARP routes return the amino acid at Hou 2026's measured 0.25–0.40 per ARP lost.</li><li>Acid yield and pKa stay at B60.</li></ul>Score Meynier (drift statement only), Hou Fig. 3, Zhou 2023 Fig. 2 (ARP alone), Cao Fig. 4c | 3.0 | Sohn 1995; Hou 2026 Fig. 3A / 3H; Zhou 2023 Fig. 2; Cao 2017 Fig. 4c; the Meynier pH statement (all SEEN; Meynier's volatile tables are NOT read: M1b item 7) | cooled pH · Meynier ±0.2 per label, Hou RMS ≤ 0.3, Zhou RMS ≤ 0.3 (CF-084 pass rules) · the three B60 pH anchors (Zhou 2023 ARP + Cys) stay inside their tolerance |
| CK-13 | 06 (b) | Engine-free: from Hou 2026's measured ARP and pH, is the acid lagged behind ARP loss beyond read error (acid per ARP lost vs time)? | 1.0 | Hou 2026 (SEEN) | acid-equivalent ratio · read error · none |
| CK-14 | 06 (c) | **done**: cluster 3 ran ARP pKa 7.07 / 7.6 / 8.18 (RMS 0.69 / 0.45 / 0.84); it fails alone. Rerun only combined with CK-12 (+0.25 h) | 0.25 | as CK-12 | as CK-12 |
| CK-15 | 08 (b) | With the CK-12 state: the engine's NF pH 5 / pH 3 ratio in the Hofmann 1998 ladder | 1.0 | Hofmann 1998 T2 / T5 (in-sample) | NF ratio · Hofmann's printed NF order (sign) · the Hofmann B2.3 rows |
| CK-16 | 22 (a) | Cluster-3 arm B (inline pH, `k_aa_frag` off) on the Davidek 2008 plateau | 1.5 | Davidek 2008 — **owned by M1b** (D6) | AA mol% · Davidek's SIDA SD · the Martins 2005 AA rows |

### Batch G — Strecker: anoxic flux, O2 dose, donor sinks (G6, G9, G11)

Run after CK-01. Coordinate with the M3 owner: the `m3-phase1-build` prereg v2.1 is a DRAFT. Declare these checks as
an amendment-listed diagnostic before the freeze; they read only SEEN Table 2 cells. Total ≈ 12.5 h work, ~1 h cpu.

| id | MF / cand | what to compute | work | spends | unit · noise band · in-sample rows |
|---|---|---|---|---|---|
| CK-50 | 33 (a), 29 (d) | Arithmetic: µmol O2 available in the vial (fill and volume from CK-01) for air vs O2-flushed, against R2's 250–700 and R1's 4–8 µmol | 0.5 | Table 2 Methods (SEEN) | µmol O2 · the vial volume's print rounding · none. Under air, R2 is refuted if available O2 < 250 µmol |
| CK-43 | 29 (a), (b), (e) | From the existing V-time draws: decompose the argon total by exit (e1 / e3 / eF) and re-evaluate with <ul><li>(a) k_d at the fed 0.008 /min (CF-121);</li><li>(b) s_DG1 at Granvogl's lower end;</li><li>(e) ARP 90 mM.</li></ul>One at a time, then together | 2.5 + 0.5 cpu | Table 2 argon (SEEN); Fig. 1 (SEEN); Granvogl (SEEN) | PA + PAA mM · ±0.07 dex (read error ⊕ 7–16 % run-to-run) · Table 2 air and air + Cu (identities) and the Fig. 1 U1e / U1l rows stay in their bands |
| CK-44 | 29 (d) | Re-solve σ_ox with the air arm as O2-flushed (×4.8 O2), then recompute k_anox and the argon total | 1.0 | as CK-43 | as CK-43 |
| CK-46 | 25 (a), 24 (a), 29 (c) | Add MGO + AA → adduct with Piloty's pH-flat consumption (≥ 3.5e-5 L/(mmol·min), a lower bound) and Kanzler's amine order 0.75–0.84. Forward on <ul><li>Piloty glycine loss;</li><li>Kanzler MGO loss with / without amine;</li><li>the M3 argon donor budget.</li></ul>No refit | 3.0 + 0.5 cpu | Piloty 1979 T1 / T4; Kanzler 2026 Figs. 1–2 (SEEN) | % Gly lost · Piloty's analyser CV (unprinted: declare ±5 %) · trunk MGO rows of B1 and B18's 10 Strecker / pyrazine rows stay in their σ_log |
| CK-47 | 24 (b) | Arithmetic: Kanzler's k2 at trunk concentrations (≤ 0.1 mM MGO); is a second-order self-loss negligible there, so that `k_mgo_mel` must be amine chemistry? | 1.0 | Kanzler 2026 (SEEN) | MGO k at dilution · Kanzler's two-point spread (×3.2 within lab) · none |
| CK-45 | 30 (a) | Forward B18's glucose-pot PA with `k_mgo_ak` at Hofmann's fed MGO + Phe measured rate in place of the fitted one (side effect on pyrazines reported) | 2.0 | hofmann2000 Tables 3–4 (SEEN, exact calibrations); hofmann2000 Fig. 2 (SEEN) | PA mM per time · the SIDA read error ±0.05 dex · B18's 10 fit rows (pyrazines) |
| CK-49 | 32 (b) | With the CK-12 pH machinery: the in-situ pH of the 0.5 M ARP trapped pot under 0.5 mol/L vs 0.5 mmol/L phosphate, and the decay slowdown it implies | 1.5 | the hofmann2000b trapped-pot text (SEEN as dossier text) | k_d /min · ±0.005 /min (10 % CV) · Fig. 1 rows. NOT RUNNABLE beyond the text: Figs. 2–3 are BLIND |
| CK-48 | 31 (b) | After CK-08: the Pan rows with a methionine Strecker step at a cited constant | 1.0 | Pan 2025 (SEEN). The Deng half is NOT RUNNABLE | aldehyde ratio · Pan's SD · B22 rows |

### Batch C — the thiol sink and the formation barriers (G3, G7)

One forward set with measured barriers serves MF-01, -10, -11 and -16. Total ≈ 9.5 h work, ~1 h cpu.

| id | MF / cand | what to compute | work | spends | unit · noise band · in-sample rows |
|---|---|---|---|---|---|
| CK-23 | 01 (a) | Bound `k_mft_decay` from US 10,076,126's co-heated labelled-standard survival (as quoted in CF-013). Forward Schieberle T4 ratios with the sink at that bound and supply as shipped | 2.5 | the US 10,076,126 standard recovery (SEEN via CF-013). Its MFT **shape** is a declared B59b hold-out and is NOT read | MFT 360/30 and 720/30 · the rows' σ_log · the B60 Hofmann 145 C rows |
| CK-26 | 11 (a); side 10 | **Owned by M2 (coordinator, 2026-11-06, from M0 batch 2, PR #89):** the Zhai 2023 JAFC T-order includes the 140 C pots (M2's T-140 / T-140c) and pot D is M2's Rt-D; without them only B60 in-sample cells remain. Forward with per-route measured class barriers (Cys → H2S 123–135, furfural 147, NF formation 81–93; atlas F23-d) instead of the single 64 kJ/mol. Score the Zhai 2023 JAFC four-time T-order (16 comparisons) and the Cerny pot D share | 3.0 + 0.5 cpu | Zhai 2023 JAFC Fig. 3 grid (SEEN); Cerny pot D (SEEN; its decisive use is **owned by M2**, so it is reported as a side reading only) | T-order sign · the grid's read error (semi-quantitative: sign only) · the B60 145 C rows stay in σ |
| CK-27 | 11 (b), 01 | **Owned by M2 (coordinator, 2026-11-06, from M0 batch 2, PR #89):** it is M2's own requirement that the T verdict hold at the Gigl-ceiling profile. As CK-26, with `Ea_decay_thiol_sink` at the median of Gigl's measured thiol-loss barriers (range 7–102 kJ/mol; read the values from the source) instead of the 102 ceiling | 1.0 | as CK-26 | as CK-26 |
| CK-28 | 16 (b) | **Held (coordinator, 2026-11-06, from M0 batch 2, PR #89):** its shares (Zhou 2023, Zhang 2024) are M1's D3 / D4 and M1's part C was never run. Engine-only: the maximum disulfide share reachable with `k_dimer_mft` at its instant bound, against the measured shares. Is the formation rate or the oxidant supply the limit? | 1.0 | Zhou 2023, Zhang 2024 (SEEN) | share · the rows' σ · none |
| CK-24 | 01 (b) | **done** (EX-061 not separated) | 0 | — | — |
| CK-25 | 01 (c) | Guth 1995 release rate in C3 — **owned by M1** (C2 / C3) | — | Guth 1995 | — |
| CK-29 | 18 (a) | Engine-free: within-study MFT drop vs printed lipid-carbonyl load (Whitfield 1988, Elmore 2002, Yiltirak 2026 arms). Is a first-order capture in the carbonyl load consistent across them? | 2.0 | the three sources (SEEN). Xu 2011 excluded (E12 TEST); E7 TEST F30#11 / 15 / 16 NOT RUNNABLE | ΔMFT dex · each study's SD · none |

### Batch E — headspace, sulfide partition, furfural lifetime (G9, G10)

Total ≈ 7 h work, ~0.5 h cpu.

| id | MF / cand | what to compute | work | spends | unit · noise band · in-sample rows |
|---|---|---|---|---|---|
| CK-34 | 13 (a) | Forward the B2.3 fed-sulfide rows with an H2S vessel partition (Kaw(T) from `henry_constants.yml`; Vg/Vl per pot, a declared band when unstated). Report the bias on `k_fur_fft_hs` / `k_nf_mft` and the Munch cysteine / H2S contrast | 2.0 + 0.5 cpu | Hofmann 1998 T3 / T4 / T10; Munch 1997 T5 (SEEN) | contrast ratio · σ ≥ 0.1 · the B2.3 fed-sulfide rows |
| CK-35 | 13 (b) | Search the dossiers for a measured furfural thermal loss in phosphate (CF-126's Bork 2024 recovery is a lead). If found, forward with `k_fur_decay` at it | 2.0 | that source if found | as CK-34 |
| CK-36 | 13 (c) | No seen data (direct FFT stability is CF-065 decider iii) | — | — | NOT RUNNABLE on seen data |
| CK-37 | 14 (a) | Engine-free: extend CF-068's O2 : (Cys + pentose) residual regression to every isotope-quantified pot with a printed headspace (US 8,920,862 level; Mottram 2002 if quantified) | 1.5 | Yiltirak 2026 (a hold-out, "not unseen"); US 8,920,862 (a B59b hold-out level: **owned by B59b**, so exclude it unless B59b releases it) | log level gap · τ_ABS-scale · none. The WUR blind samples are NOT RUNNABLE |
| CK-38 | 14 (b) | No seen data; an author question (aliquot and work-up) | — | — | — |
| CK-39 | 07 (a) | Forward B60 with `k_nf_decay` at Bork's measured 4.9e-3 /min (130 C, pH 5) on the in-sample rows, to show the side effects | 1.5 | Bork 2022 / 2023 — **owned by M1b(ii)** (carried as a band) | k /min · ±0.03 dex (two-reader, weighted fit) · the B60 NF / MFT rows |

### Batch B — the enolisation pH law (G2)

Total ≈ 4.5 h runnable work plus 5 h owned by M1b, ~0.5 h cpu. Several checks read M1b's decisive DEV data and are
owned by it.

| id | MF / cand | what to compute | work | spends | unit · noise band · in-sample rows |
|---|---|---|---|---|---|
| CK-18 | 09 (b) | Engine-only attribution: split the engine's FFT pH 4.5 → 6.5 fall into the furfural-supply factor and the FFT-step speciation factor | 1.5 | none (engine only) | — (a diagnostic that sends CK-17 to the right step) |
| CK-17 | 09 (a) | Forward with an amine-catalysed pentose enolisation law rising with pH, at Martins' slope on the mean pH (0.55 dex/pH, CF-118). Score the FFT / furfural / MFT ladders of Cerny 2007 and Xu 2007; side reading on Hofmann T2 | 3.0 + 0.5 cpu | Martins 2003 (constant source) — **owned by M1b** (D5); Cerny 2007, Xu 2007 (SEEN) | log ratio per pH pair · ω 0.42 · the Hofmann 1998 T2 rows |
| CK-19 | 19 (a) | `TRUNK_PH_TERMS` on, at the mean-pH slope; forward on the desclaux2006 / davidek2008 flags | 2.0 | Davidek 2008 — **owned by M1b** (D6). Meynier NOT READ (item 7) | the flagged rows' σ · the R3En trunk rows |
| CK-20 | 20 | No seen-data discriminator (DM entry: "none is on disk") | — | — | NOT RUNNABLE on seen data |
| CK-21 | 28 (b), (c) | After CK-06: forward the Schibilsky pots with a 2,3-enolisation pH law at a cited slope, ± HDMF phosphate ×4.3 | 2.0 | Schibilsky 2019 (SEEN) | DMHF ratio pH 8 / 5 · Schibilsky's SD · the acrylamide-lane fit rows |
| CK-22 | 17 (b) | Arithmetic: does Xu 2007's measured ACTZ pH slope (−0.26 / −0.41) match an NH3-fraction-driven slope, or the opposite? | 1.0 | Xu 2007 (SEEN) | log ratio · ω 0.42 · none |

### Batch F — the route shares and the donor dose

Total ≈ 3 h work.

| id | MF / cand | what to compute | work | spends | unit · noise band · in-sample rows |
|---|---|---|---|---|---|
| CK-42 | 15 (a) | **Done (coordinator, 2026-11-06, from M0 batch 2, PR #89):** the H2S half is refuted (0.65 mM Cys releases 2-3 uM sulfide in 2 h, >= 17x below the Table VI H2S arm). Also spends v2 TEST(E18) F27#15-17, already scored in W1.1. See CK-42b. Engine-free: Grosch Table III's optimum-then-fall in Cys dose. Is it consistent with Cys → H2S at the measured Cys rate plus a Cys-dose-dependent MFT capture (the CF-065 cysteine depression)? | 2.0 | Grosch 1992 T2 / T3 (SEEN, scored) | MFT ratio · Grosch's SIDA SD (±1.2 / 3.7 µg) · none |
| CK-42b | 15 (a) | **Done (PR #89):** the Grosch cysteine term follows total free-thiol concentration, not the cysteine : thiamine ratio (the optimum stays ~340 umol when thiamine drops to 6 umol); a total-thiol law built on Table III predicts the GSH pots unfitted (no. 7 -0.11 dex, no. 8 +0.11); the Cys.HCl pH artefact is retired (delta-pH -0.13). Nominates a thiamine-route factor g(total free thiol) for an MF-15 case; decisive: CF-089's designed series with thiamine measured plus a non-amine thiol arm | 2.0 | Grosch 1992 T2 / T3 / VI; v2 TEST(E18) F27#15-17 (already scored; now SPENT for the MF-15 case) | as CK-42 |
| CK-40 | 10 | **owned by M2**. CK-26 supplies a side reading only | 1.0 | Cerny pot D, Givaudan | — |
| CK-41 | 12 | **done** (atlas sec. 5: `HEXOSE_C5_ENTRY` gives 45–78, or 0.79 with Gly). The rest is **owned by M4** | 0 | — | — |

### Batch H — pyrazine ring step

| id | MF / cand | what to compute | work | spends | unit · noise band · in-sample rows |
|---|---|---|---|---|---|
| CK-52 | 35 (a) | Forward with `k_cond` at the analogue's measured constant (Ea 44, E14); score the Leahy T3 ratio and the CF-087 R rows | 2.0 | Leahy (B18 fit rows, in-sample, cannot decide per M6); Zhou 2023c T1 (SEEN). E13 / E14 TEST NOT RUNNABLE | log ratio · read error · B18's 10 rows |
| CK-51 | 34 | **owned by M6** (DM:pyrazine_ring_aminoketone) | — | — | — |
| CK-53 | 36 | **done** (atlas sec. 5) | 0 | — | — |

### Batch I — buffer catalysis (G5)

| id | MF / cand | what to compute | work | spends | unit · noise band · in-sample rows |
|---|---|---|---|---|---|
| CK-54 | 26 (a), 05 (a) | Which steps must carry κ_P (DEV band [5, 40] /M) so that Mottram's MFT moves ×2.2–2.6 and FFT does not? | 3.0 | Mottram 2002, Davidek 2002 — **owned by M5** (D1, D3) | — |
| CK-55 | 26, 32 | **done** (CF-120: weak support at pH 7.5 only) | 0 | — | — |
| CK-33 | 05 | **owned by M1c / M5**; no seen-data discriminator | — | — | — |

### Batch J — glycation

| id | MF / cand | what to compute | work | spends | unit · noise band · in-sample rows |
|---|---|---|---|---|---|
| CK-56 | 27 (a) | Arithmetic: Wolf 1977a's E recomputed with the plateau end points excluded from k | 1.5 | Wolf 1977a (SEEN) | E kJ/mol · Wolf's residual SD (±8.4 kJ/mol) · B20 rows. It explains if E falls within ±8.4 of 96.8 |

### Batch K — acrylamide (after CK-04)

| id | MF / cand | what to compute | work | spends | unit · noise band · in-sample rows |
|---|---|---|---|---|---|
| CK-57 | 37 (a), (b) | Forward Chang 2021 with its printed heat-up (if CK-04 finds one), and with the phosphate pH law off for its water arm | 3.0 | Chang 2021 (SEEN) | acrylamide ratio 30 / 10 min · Chang's SD · the B3 fit rows |
| CK-58 | 38 (a) | A fructose + Asn step at the measured Fru / Glc ratio (×1.1–2.9, F29-a, the band's lower end first); forward Lin 2022 | 1.0 | the F29-a sources (constant), Lin 2022 (SEEN) | acrylamide · Lin's SD · the B3 rows |
| CK-59 | 39 (a) | After CK-04: the precursor charge from the SPI composition (Gorissen 2018) | 1.5 | Gorissen 2018; the SPI row (SEEN) | acrylamide · the row's SD · none |
| CK-60 | 40 | No seen-data discriminator (CF-106 decider) | — | — | NOT RUNNABLE on seen data |
| CK-61 | 41 | **done** (CF-104 T4 – T6) | 0 | — | — |

### Batch L — lipid (after CK-05)

| id | MF / cand | what to compute | work | spends | unit · noise band · in-sample rows |
|---|---|---|---|---|---|
| CK-62 | 42 (a), (b) | Forward the three pots with the matrix's printed hydroperoxide / PUFA charge (where printed) vs the default | 2.0 | li_2026, bi_2020, Trikusuma (SEEN) | log level · τ_lipid (lane) · the lipid-lane fit rows |

### Not runnable, and checks owned by a case

| reason | checks | dataset |
|---|---|---|
| BLIND | CK-49 beyond the text; any decisive use for MF-29 / MF-32 | hofmann2000b Figs. 2–3 (M3 V-blind) |
| BLIND | CK-37 (WUR part) | the WUR blind-test samples |
| reserved / TEST | CK-08 (Deng half), CK-48 (Deng half) | Deng 2022 (M3 phase 2) |
| reserved / TEST | CK-29 (E7 part) | E7 TEST F30#11, #15, #16; Xu 2011 (E12 TEST) |
| reserved / TEST | CK-51, CK-52 (TEST part) | E13 / E14 TEST studies (owed, unscored) |
| reserved / hold-out | CK-23 (shape), CK-37 (US 8,920,862) | US 10,076,126 shape; US 8,920,862 level (B59b hold-outs) |
| no seen-data discriminator | CK-20, CK-36, CK-60 | — |
| owned by M1c | CK-30, CK-31, CK-32, CK-33 | Zhai 2021 Fig. 4b, Zhang 2025b blank, Cao 2017, Zhao 2017b |
| owned by M1b | CK-16, CK-17, CK-19, CK-39 | Davidek 2006a / 2008, Martins 2003 / 2005b, Bork 2022 / 2023, Mikami / Nakamura |
| owned by M1 / M2 / M4 / M5 / M6 | CK-25, CK-40, CK-41, CK-54, CK-51 | Guth 1995; Cerny pot D; Hofmann glucose pots; Mottram / Davidek 2002; aminoketone pots |

### Total hours

Runnable now (not owned, not done), single core:
- batch 0: 9.25 h;
- batch A: 5.25 h (without CK-16);
- batch G: 12.5 h;
- batch C: 9.5 h;
- batch E: 5.5 h (without CK-39; CK-37 partial);
- batch B: 4.5 h (CK-18, -21, -22);
- batch F: 2 h;
- batch H: 2 h;
- batch J: 1.5 h;
- batch K: 5.5 h;
- batch L: 2 h.

That is **≈ 59.5 h of work plus ≈ 3 h of unattended compute.** Write-ups add about 15 %, so call it **70 h (range
55–90 h)**.

The case-owned checks add ≈ 12 h inside their own preregs (CK-16 1.5, CK-17 3, CK-19 2, CK-39 1.5, CK-54 3,
CK-40 1).

Batch 0 should run first. Historically integrity is the most common cause here, so batch 0 is likely to retire some
misfits outright and shrink the later batches.

---

## 4. Shared-structure groups

Each group is a set of ≥ 2 misfits that one shared piece could explain. P is my rough probability that the named
piece explains at least two members once checked; it is not a probability that the piece exists.

| group | piece | members | P | the batch that tests it once |
|---|---|---|---|---|
| G1 | the in-situ pH state: amine fate and acid timing in `ph_state` | MF-06, -08, -22, -32 (acid load); side: -03, -17 | 0.6 | A (CK-12) |
| G2 | the enolisation pH law (`acid` = [H+] on the pentose / ARP steps; pH-free trunk steps) | MF-09, -19, -20, -28, -36; side: -08 | 0.5 | B (CK-18 → CK-17) |
| G3 | the thiol-sink law: frozen, pot-independent `k_mft_decay`, Ea on its ceiling | MF-01, -11, -16, -18; side: -14 | 0.5 | C (CK-23, CK-27) |
| G4 | the supply reservoir (TTCA / CAMA; instant supply, CF-107 compensation) | MF-01, -02, -03, -04, -05, -15 | 0.55 | M1c (owned) |
| G5 | a step-specific buffer factor (phosphate κ_P) | MF-05, -26, -32 (CF-120), -37 (buffer law), -28 | 0.45 | I (owned by M5) |
| G6 | dicarbonyl–amine adducts (a missing donor and amine sink) | MF-24, -25, -29; sign check on -34 | 0.35 | G (CK-46) |
| G7 | per-route formation barriers (E11) | MF-10, -11 | 0.6 | C (CK-26) |
| G8 | dicarbonyl / osone supply from hexoses | MF-12, -23, -30, -34 | 0.3 | G (CK-45) plus M4 / M6 |
| G9 | a vessel and headspace operator (O2 dose, Vg/Vl, H2S partition) | MF-13, -14, -29, -33; side: -07 (O2 / metal), -23 | 0.4 | G (CK-50, CK-44) + E (CK-34, CK-37) |
| G10 | furfural lifetime and sulfide speciation in the FFT step | MF-09, -13 | 0.35 | B (CK-18) + E (CK-34, CK-35) |
| G11 | the Amadori (ARP) decay law (order, vessel, acid) | MF-32, -29, -06 (acid downstream of the ARP) | 0.35 | G (CK-43, CK-49) + A (CK-13) |
| G12 | integrity: a keying, unit or pot-definition error | MF-21, -23, -28, -31 (Pan), -39, and others never re-read | 0.5 that ≥ 2 of the 16 never re-read retire at step 1 | 0 |

Most misfits fall in more than one group. Because they overlap, run batch A before batch B: an enolisation pH law
must be judged on the in-situ pH, not on the label. Likewise run batch G's O2 arithmetic (CK-50) before any M3 rival
is credited with an O2 route.

---

## 5. The top checks, ranked

Score = (leverage × P(discriminating)) / cost.
- **Leverage** counts the W5.4 cells and the misfits the check touches.
- **P(discriminating)** is my estimate that the check separates its candidates, or retires the misfit.

| rank | check | what | cost (h) | data spent | why |
|---|---|---|---|---|---|
| 1 | CK-06 | Schibilsky: is the pH passed to the lane at all? | 0.75 | Schibilsky (SEEN) | Identical outputs at two pH values is the pattern of an input bug. P 0.5; it retires MF-28 or confirms a structural pH gap |
| 2 | CK-01 + CK-50 | Re-read Hofmann 2000b Table 2 and Methods; O2 arithmetic | 1.25 | Table 2 (SEEN) | It touches MF-29, MF-33 and the U3a denominator of M3's primary claim. P 0.45. Cheap, and it must precede the M3 freeze |
| 3 | CK-12 | pH state with the amine returned (Sohn, Hou), forward | 3.0 | Sohn, Hou, Zhou 2023 Fig. 2, Cao Fig. 4c, Meynier pH statement (SEEN) | G1 reaches the pH lever for every odorant and four misfits. P 0.55 |
| 4 | CK-02 | zhou2025b −16σ rows | 0.75 | zhou2025b (SEEN) | A floor or zero prediction is likely (P 0.6); it cleans the R3En trunk flags |
| 5 | CK-43 (+CK-44) | Argon flux decomposition by exit; fed k_d, s_DG1, ARP purity, O2-flushed air | 3.5 | Table 2, Fig. 1, Granvogl (SEEN) | It names the candidate MF-29 lacks; G11 and G9. P 0.4 |
| 6 | CK-46 | MGO–amine adduct sink at Piloty's / Kanzler's measured rates | 3.5 | Piloty, Kanzler (SEEN) | G6: MF-24, -25, -29; Strecker and pyrazine levers × cys / sug. P 0.35 |
| 7 | CK-26 (+CK-27) | Per-route measured barriers; sink Ea at the median of Gigl's measured barriers | 4.5 | Zhai JAFC grid, Cerny pot D (side reading) | G7 + G3: MF-10, -11, -01; Th × Tt, Th × ck. P 0.4 |
| 8 | CK-23 | Sink rate bounded by the co-heated standards; Schieberle ratios | 2.5 | US 10,076,126 standard recovery (SEEN); not its shape | MF-01 candidate (a), against M1's open non-separation. P 0.35 |
| 9 | CK-34 | H2S vessel partition on the fed-sulfide rows | 2.5 | Hofmann T3 / T4 / T10, Munch T5 (SEEN) | G9 / G10: MF-13 and the bias on `k_fur_fft_hs`. P 0.4 |
| 10 | CK-42 | Grosch Cys-dose optimum, engine-free | 2.0 | Grosch T2 / T3 (SEEN) | MF-15: the Th × thi and Th × cys levers. P 0.35 |

---

## 6. Misfits where step 1 (integrity, R) was never done

Integrity has been the most common cause in this repository, so these misfits come first. Each has keyed rows (K) at
most, and no re-read because of the misfit.

| MF | what | check | what to read |
|---|---|---|---|
| MF-29 | argon PA + PAA ~5x high | CK-01 | Table 2 basis, vial headspace, deoxygenation, air vs O2, the phosphate print |
| MF-33 | O2 dependence (R2 O2 demand) | CK-01, CK-50 | the same Methods paragraph |
| MF-31 | Pan 2025 −3.6 to −5.6 dex | CK-08 | units and basis (Deng half reserved) |
| MF-21 | zhou2025b xylose-alone −16σ | CK-02 | keyed rows and the engine output |
| MF-23 | Steinhagen osones | CK-03 | autoclave profile, O2, units |
| MF-28 | Schibilsky furanics, identical outputs | CK-06 | whether the bundle passes pH |
| MF-37 | Chang 2021 / Ye 2024 acrylamide | CK-04 | heat-up, basis |
| MF-38 | Lin 2022 fructose | CK-04 | units, sugar charge |
| MF-39 | extruded SPI −3.98 dex | CK-04 | units, precursor charge |
| MF-42 | lipid matrix rows | CK-05 | extraction, headspace, units |
| MF-17 | ACTZ −3.25 dex | CK-07 | hold-out status (Amendment 45 vs the panel); CF-009 partial |
| MF-16 | disulfide share | CK-09 | the species and basis of the share |
| MF-15 | Grosch Cys dose | CK-10 | the engine's pot charges |
| MF-10 | Givaudan ribose / NF ratio | (CK-40, owned by M2) | the Givaudan cells |
| MF-34 | pyrazine total (Leahy T4) | CK-08 | Leahy rows (CF-077 fixed the code, not the source) |
| MF-30 | glucose-pot PA lag | — | hofmann2000 Fig. 2 (EX-049 keyed; no re-read for the lag) |

---

## 7. How to keep this page

- **One writer at a time.** New misfits take the next MF id; check against `origin/main` and open branches first.
- **A misfit leaves this page** only when its CF entry (or EX entry) records a resolved status. Move it to a "closed"
  list with that id; do not delete it.
- **Record each run in the checks table.** A run check records its outcome (explains / refutes / inconclusive) and
  its record path in its row. Its *spends* column then becomes binding on every case.
- **When a case or prereg freezes,** compare its decisive data with the *spends* column. Any overlap goes into the
  case's SEEN list.
- **Open items found while compiling this page** are left to their owners:
  - **No CF entry exists for MF-29**, the argon over-prediction; propose one.
  - DM:strecker_routes does not mention MF-29.
  - The Zhang 2025b ACTZ hold-out status (CK-07).
  - The B60 profile file's header still reads "B8 sulfur objective".
