# Case M1c: the thiol supply structure (a base-favoured TTCA exit and Zhang 2025b's missing pool)

*Plan step M1c ([STATUS_AND_PLAN sec. 5, W2](../../guides/STATUS_AND_PLAN.md), as re-sequenced 2026-11-05; error E10,
with E12 and CF-076 / CF-063 / CF-107). Items 1 - 5 of sec. 7, plus the bands every new coordinate carries, written
2026-11-05 before any result of the change exists and before any newly downloaded thiol-supply paper was read. No
engine was run for this case, on any pot. Items 6 - 9 are appended later and never edited into this part.*

**Blinding.** Eleven newly downloaded PDFs bear on thiol supply: `canle1996`, `gajda1990`, `belikov1967`,
`wlodek1993`, `davidek2003`, `davidek2004`, `carta1996`, `szilagyi1979`, `davis1991`, `davis1991b`, `sohn1995`. None was
opened, and no extraction of them was read. What this case says about them comes from the private download queue only
(titles; abstract-level descriptions for some, title only for `belikov1967`, `carta1996`, `szilagyi1979` and
`sohn1995`). Sec. 4.3 fixes now which of them could be decisive and on what measurement, and which are band sources.
A paper named as decisive can never set a band, and a band source can never decide.

**Seen by the author.**
- **Engine results:** the supply screen's outcomes as the parked M1 + M2 prereg v2 records them (sec. 0B table:
  every structure's feasibility on signatures (a) - (c), including engine predictions on Zhang 2025b's blank pot:
  best TTCA 11.6 / 8.8 mmol/L under S1, free Cys >= 3.2 / ~10 / 15 mmol/L when TTCA and CAMA pass under S2 / Kallen /
  S3); the engine's cooled-scale pH in Cao 2017's pot (~0.2 low at 30 min under C1); M1's decisive comparison and
  results review (EX-061: R0 and C1 per-species residuals on Zhang 2025b and Schieberle 2000); the no-refit probe
  numbers on Zhai 2019 Fig. 3A (TTCA -0.29 dex mean, CAMA +0.71 dex) and B59b's prereg check (CAMA formation 2 - 5x
  too fast), relayed from the repository's records by a reading agent.
- **Measured values only (no engine number):** the dossiers named in sec. 2, the F04 atlas file's summary lines, and
  the queue's 2026-11-05 section.
- **Not seen:** any engine result on Hou 2017, Wei 2022, Hou 2026, Zhai 2020, Kang 2026's storage and pH-ladder
  tables, or Mizutani 1960 (E10's only v2 TEST contrast, F04#26: read in item 7 only). The screen's own record
  (`kinetic_core_supply_structure_signature_screen.md`) did not exist when this was written; its numbers were taken
  from the prereg's sec. 0B and the plan row.

## 1. The claim

In a heated cysteine + pentose pot, cysteine is held in a reservoir (TTCA and the cysteine Amadori compound, CAMA), as
in M1's C1. The claim adds two things. First, TTCA leaves the reservoir through an **irreversible exit that is faster
at higher pH**, so that fed TTCA is lost ~10x faster at label pH 8 than at 5.5 (Zhai 2021 Fig. 4b). Second, Zhang
2025b's ~19 mmol/L of each reactant that no measured pool holds is **accounted for by named pool members**, each
entering with a measured or same-class constant or a declared band, never as one free pool sized to the gap. The size
of the account is an output. If the members' bands cannot reach it, the claim says so rather than widening them.

Network steps that change (names from [sulfur.py](../../../src/kinetic_core/sulfur.py) and
[parameters_sulfur.py](../../../src/kinetic_core/parameters_sulfur.py)):

| part | step(s) | shipped (R0, B60) | C1 (M1) | M1c |
|---|---|---|---|---|
| A. the reservoir | `r_ttca_form` / `r_ttca_cys` (Cao's K, detailed balance); `r_cama_form`; `r_cama_dpo` / `r_cama_tdp` | off (`TTCA_FORMATION`, `CYS_AMADORI_ROUTE` disabled): instant supply through `r_pent_dpo` / `r_pent_tdp` | on; K 0.060 L/mmol, pH- and T-flat; `r_cama_form` pH-flat (3.29e-3 /min at 90 C, anchored on Cao's blank, Ea 100.9 kJ/mol); CAMA decay on de Roos 2005 Fig. 4 | as C1, with the release `k_ttca_cys` free (M1 review: C1's reservoir is 3 - 5x too slow) |
| B. the base-favoured exit | a factor 10^(n (pH_c - 5.5)) on ONE irreversible TTCA exit. Variant **X-cama**: on `r_cama_form`. Variant **X-deg**: on `r_ttca_deg` | none | none | both variants built and both run (no choice between them on scores). n on the cooled-scale pH (sec. 5.2) |
| C. pool member P1: open-chain / Schiff-base TTCA | a new species TTCA_oc in fast ring-chain exchange with TTCA, invisible to Zhang's MRM and ELSD; it reacts as TTCA does | absent | absent | f_oc(T) = [TTCA_oc] / [TTCA] inside its band |
| D. pool member P2: the minor C-2 epimer | observation operator only: a share of TTCA that Zhang's MRM may not integrate | absent | absent | share inside its band (a declared operator switch, no new chemistry) |
| E. pool member P3: cystine | Cys -> 1/2 cystine, first order, a new step (not B11's `ch_cys_ox`, which needs OX) | absent | absent | a ceiling rate from Hou 2026. It takes cysteine only, so it cannot by itself make a 1 : 1 deficit (CF-063) |
| F. pool member P4: the cysteine Amadori compound onward | a second CAMA isomer outside Zhang's MRM transition; the 3-deoxypentosone + cysteine adduct `r_tdp_cys_trap` (1 Cys : 1 C5) | absent / zero | absent / zero | numerical bands only (sec. 5.2): a verdict that needs P4 is "conditional on a numerical band" and does not ship |
| not changed | Cao's K level; the CAMA decay law; the thiol loss law (`k_mft_decay` as M1 left it); the pH balance (M1b(i)) | | | |

The run uses the measured pH trajectory where one exists (Cao 2017; Kang 2026's end points for the Zhai-design pots,
sec. 5.2), whatever M1b(i) finds. Every pH comparison states its scale (cooled or in situ).

## 2. All the evidence, by kind

"Engine": whether the measurement has been run under the engine (**seen**) or not (**not run**). Laboratories as in
[lab_ids.yml](../contrast_benchmark/lab_ids.yml). For E10 the frozen v2 split
([split_v2_per_change.yml](../contrast_benchmark/split_v2_per_change.yml)) has DEV labs BTBU (Xie), the Ho network
(Jiangnan) and Kallen (Pennsylvania), and one TEST contrast, F04#26 (Mizutani 1960, Hokkaido): item 7 only.

### 2.1 Isotope labelling

| source | what it shows | for / against | engine |
|---|---|---|---|
| Zhai 2021 Fig. 7a (13C5-xylose + fed TTCA; CF-091) | label crosses into TTCA: release and recombination occur at 100 C. Quantified only as a difference of two pots; the SI shows the m/z 259 ion but not its share | for a reversible release (C1, X); at face value it disfavours CF-076's E4c point (z 6.3) | seen (CF-076 cluster 5) |
| Zhai 2023 Food Chem Table 2 (13C5-xylose + fed TTCA, 120 min) | exogenous xylose carbon in thiophenes 0 % (100 C), 19 - 21 % (120 C), 34 - 45 % (140 C); no n | release + recapture grows with T | not run (thiophenes are not scored) |
| Zhao 2019 (BTBU; [13C5]-TTCA, 140 C, pH 5.5, 30 min) | 97 % of MFT carries the intact TTCA C5 when TTCA dominates the pool | TTCA's C5 does reach MFT (C1, X; R0 has no TTCA in such pots) | not run |
| none | no labelling of the base-favoured loss's products, and no label on the missing pool | gap | -- |

### 2.2 Fed intermediates

| source | what it shows | for / against | engine |
|---|---|---|---|
| **Zhai 2021 Fig. 4b** (Jiangnan; TTCA 10 mM +- Xyl, unbuffered, label pH 5.5 / 6 / 7 / 8, 100 C) | TTCA left at 80 min 0.819 (5.5) vs 0.147 (8): ln-fraction lost 9.6x; 0 - 40 min first-order equivalent 0.0017 -> 0.0164 /min (~0.4 dex / pH). Every curve ACCELERATES (pH 7: 0.0045 /min over 0 - 20 min, 0.049 over 80 - 100). Error bars only on pH 5.5 / 6 / 7 TTCA alone | for part B. The acceleration is **against a plain first-order exit on a pH that drifts acid** (it would decelerate; see Kang below) and for an induction: an open-chain pool that builds up (P1), or autocatalysis | **seen** (signature (b); B16 / route rows; CF-076) |
| Zhai 2021 Fig. 5 (same pots) | standing 3-DX per TTCA lost at 120 min: 19 / 17 / 5 / 0.8 % (pH 5.5 / 6 / 7 / 8); 1-DX 32 - 106 umol/L; A420 doubles 5.5 -> 8 | the fast base-favoured loss goes **neither to 3-DX nor (measurably) to 1-DX**: into browning or an unmeasured product. Bears on where X-cama / X-deg send the carbon | seen (deoxyosone sub-model rows) |
| Zhai 2021 Fig. 3a (pH 5.5, 100 / 120 / 140 C) | zero-order k 0.0274 / 0.0651 / 0.0813 mmol/L/min; printed Ea 80.99 not reproduced (~35 kJ/mol) | the exit's temperature law at one pH | seen |
| Kang 2026 main text p. 3242 (Jiangnan; TTCA 10 mM, the Zhai design) | **measured** post-reaction pH 5.1 / 4.9 / 4.7 from label 5.5 / 7 / 8 (cooled); temperature of this series unresolved (100 or 120 C) | the Zhai-design pots end acid, so the label-pH contrast is a contrast in the first minutes only. With Zhai 2021's acceleration, a first-order base-favoured exit on this trajectory has to explain a loss that speeds up while the pH falls | not run (an input, sec. 5.2) |
| Kang 2026 SI Fig. S3 / S4 | free Cys released from fed TTCA, peak <= 16.3 mol % at 100 - 140 C; the F04 atlas records Fig. S3 as the same data set as Zhai 2021 Fig. 6 (not independent) | for a release that recombines (C1) | seen (B2.1 fits) |
| Kang 2026 SI Table S5 (pH ladder, FFT / MFT 2.91 / 3.02 / 2.73) | FFT / MFT pH-flat over label 5.5 - 8 | weak: S5's pH-7 column is S4's 100 C column (dossier sec. 4c) | not run; compromised |
| Hou 2017 (BTBU; fed TTCA, fed CAMA, fed Gly-ARP, 30 mM, 0.2 M phosphate pH 5.5, 140 C, 5 - 30 min, pH and UV followed) | reactivity Gly-ARP > CAMA > TTCA; both TTCA alone and CAMA alone make MFT and FFT; volatiles as peak areas | for TTCA's thiols passing through a slower step than CAMA's (C1, X-cama); one pH | not run; **unkeyed** (the dossier carries no numbers) |
| de Roos 2005 Fig. 4 (Givaudan; fed CAMA, boiling 0.5 M phosphate, pH 3 - 6) | CAMA decay U-shaped, slowest at pH 5 (4.5e-3 /min), 1.43e-2 at pH 6 | already the engine's CAMA decay law | spent (B59b constant; keyed held-out for E10, never a fit row) |
| de Roos 2005 text | "At pH < 7, the route to Amadori compounds is preferred"; above pH 9 the equilibria shift toward the thiazolidine | **against X-cama** (an Amadori exit favoured by base); qualitative, no numbers | not run |
| Zhai 2020 (Jiangnan; spray-dried TTCA intermediate) | TTCA -> CAMA conversion favoured at pH 7.5 - 9.5 (CAMA 6 % -> 47 % at pH 9.5, 190 C inlet) | **for X-cama**, in opposition to de Roos's sentence; but low water and dehydration, outside the regime | not run |
| Zhai 2019 Fig. 3B (fed TTCA 5 mM +- Cys, 120 C) | dicarbonyls only; no TTCA time course | context | seen |
| Zhai 2023c (TTCA vs Cys-ARP 50 mM, pH 7.4, 120 C) | no measured TTCA, ARP, Cys or Xyl time course in the main text; "complete consumption at 130 min" cited to a paper that lacks it (UNVERIFIED) | none usable | not run |
| Wei 2022 (Jiangnan; GX-ARP 5 + Cys 5 mM, water, label pH 7, 110 / 120 / 130 C) | at 120 C, 100 min: 4.24 mM Cys consumed, 0.48 mM TTCA formed (11 %) | **most consumed cysteine is in no measured pool where little free xylose exists**: for a cysteine-only sink (R-sink's cysteine arm), against a 1 : 1 adduct as the only missing form. Confounded by the GX-ARP chemistry (M1b) and an unmeasured pH | not run; unkeyed |

### 2.3 Trapping

| source | what it shows | for / against | engine |
|---|---|---|---|
| Zhai 2019 Flag 1 (120 C, TTCA + Cys) | "TTCA disappeared completely" within 10 min and a xylose peak appeared | against plain mass action; one reading is thiol addition to the open-chain Schiff base or trans-thiazolidination (a reactive open-chain member, P1) | not run; qualitative |
| none | no in-situ trap of free cysteine, of an open-chain form or of cystine at process temperature | gap: the quench itself may move the exchange (CF-076's E6, rejected on its dose band) | -- |

### 2.4 Measured intermediates and lifetimes

| source | what it shows | for / against | engine |
|---|---|---|---|
| **Cao 2017** (BTBU; Cys + Xyl 0.4 M each, 0.2 M phosphate, 90 C, 0.5 - 5 h) | K_app = [TTCA] / ([Cys][Xyl]) 58.6 / 64.6 / 60.2 / 59.0 / 58.5 / 61.1 L/mol while the cooled pH drifts 4.67 -> 4.15; pool over-closes 6 - 17 %; CAMA 0.006 -> 0.133 mol/L | a reservoir in pre-equilibrium (for C1, X); K flat under drift (against R-K; the slope of log K on pH is +0.013, 90 % -0.155 to +0.180 on measurement noise alone) | **seen** (signature (a); B59b anchor; M1) |
| **Zhang 2025b** (BTBU; Cys + Xyl 40 + 40 mM, 0.2 M phosphate pH 6.5, 120 C, 0 - 50 min) | TTCA plateau 17.7 / 18.2 (30 / 50 min), CAMA 1.61, free Cys 0.80 and Xyl 1.22 at 50 min: ~19 mmol/L of each reactant in no measured pool, the deficit equal for Cys and Xyl at every time and ~12 - 14 mM already at 5 min (CF-063). Cys and Xyl print one identical ELSD curve (y = 263.28 x - 58.733) | for "a 1 : 1 pool or an analytical cause"; against cystine or H2S alone (cysteine only) and sugar loss alone (xylose only) | **seen** (signature (c); M1's reservoir block; the screen) |
| Zhai 2019 Fig. 3A (Jiangnan; Xyl 165 + Cys 82.6 mM, water, label pH 7.4, 70 - 130 C) | TTCA 14.5 mM before heating; TTCA drains above ~110 C (130 C: 45.7 -> 6.9 mM); CAMA <= 5.6 mM; NMR finds no open-chain Cys-Amadori | the reservoir's T-law; no free Cys, no pH | seen (route rows; B59b; CF-076) |
| Zhai 2023 Food Chem Fig. 3a; Zhai 2023 JAFC Fig. 5 | regenerative free Cys from fed TTCA plateaus at ~1.0 mM (~10 % of the charge); fed-TTCA CAMA peaks 1.93 mM (120 C, 100 min) | release small, CAMA small in water | seen (route rows) |
| Radomski & Temeriusz 1989 (Warsaw; 15 sugar thiazolidines incl. TTCA, water, 37 C) | TTCA 2R : 2S = 4.06 at equilibrium; C-2 epimerisation (so ring opening) relaxes in minutes (k ~0.2 /min, 37 C, glyceraldehyde ring); the acyclic form is not seen by NMR | P1 small at <= 37 C; P2 at most 20 % of TTCA | band source (CF-076 battery's pKa band) |
| Karpinski 1990 (Warsaw; voltammetry, deaerated, 25 C) | release of the L-arabinose (pentose) ring AThz: 1.8e-4 /min at pH 6.0, 4.8e-5 at pH 7.4 (RSD up to 30 %); no open-ring thiol detected (LOD ~1 uM) | **release is NOT base-favoured above pH 6** (x3.75 slower at 7.4): so the base-favoured loss of (b) cannot be the release (`r_ttca_cys`) itself; it must be an irreversible exit (part B) | used as a band anchor in CF-076's resolve2 battery |
| Bustillo Trueba 2018 (Leuven; Cys + furfural / benzaldehyde / methional rings, 20 C) | ring loss bell-shaped in pH, maximum 4.4 - 5.2, slower at 6.0, nearly stable at 9.0 | as Karpinski: release is not base-favoured | not run |
| Pesek & Frost 1975 (aliphatic TCAs, D2O, T not stated) | unchanged for days in neutral D2O and in 8 M HCl; decompose only in > 3 M NaOH | as above | not run |
| Fife 1991 (2-aryl thiazolidines, 30 - 90 C) | ring opening (C-S) fast and reversible; hydrolysis of the iminium rate-determining, general-base catalysed (phosphate included) | the form of a release law; nothing transfers to TTCA | not run |
| Kallen 1971 / 1971b / 1973 (Pennsylvania; Cys + formaldehyde, 25 C) | K_app rises ~20x from pH 5.6 to 7.5 (the cysteine monoanion carries the ring); k_r ~1e-6 /s, rising from pH 2.8 to 5.6; phosphate catalysis saturates (x3.5 - 5) | the shape R-K borrows. A K that rises with pH stabilises TTCA at pH 8, **against R-K on (b)** unless an exit outruns it | seen (CF-076 battery) |

### 2.5 Mass balance

| source | what it shows | for / against | engine |
|---|---|---|---|
| Zhang 2025b (above) | 52 % of each reactant accounted at 50 min | the gap to be explained | seen |
| Cao 2017 (above) | over-closes 6 - 17 %; ceiling on all cysteine-destroying channels ~4e-4 /min at 90 C, pH 4.2 - 4.6 | against a large sink at 90 C | seen |
| Hou 2026 (Jiangnan; Xyl 0.2 M + Gly 0.1 M products with 5 mM Cys, 120 C step) | after the 120 C step TTCA 0.193 mM (3.9 % of the Cys) and cystine ~0.20 mM (~8 % as Cys equivalents): most of the cysteine is in neither | cystine is a real but minor member at 120 C (P3's band); for a cysteine fate outside TTCA and cystine | not run |
| Wei 2022 (above) | 11 % of consumed Cys as TTCA | as above | not run |

### 2.6 One factor varied inside one study

| source | factor -> result | for / against | engine |
|---|---|---|---|
| Zhai 2021 Fig. 4b | label pH -> TTCA loss x10 (2.2 above) | for part B | seen |
| Zhai 2021 Fig. 4b +- Xyl; Fig. 6 | extra xylose speeds loss modestly at pH 5.5 - 7, slows it at pH 8 in the first 40 min | release + recapture (C1, X) | seen (CF-076) |
| Kang 2026 SI Table S6 (storage, 60 days) | TTCA retention 95.51 / 96.70 / 88.81 % at pH 5.5 / 7 / 9 | non-monotone, storage at <= 40 C: a storage cell never sets a processing law | not run; context only |
| Zhang 2025b lysine arms | lysine drains TTCA and CAMA; high-lysine free Cys plots below the ELSD curve's zero-signal point 0.223 mmol/L | for a low-end ELSD problem (R-cal) | lysine is not in the sulfur lane: not representable |
| Cao 2017 glycine arms | glycine speeds the TTCA drain | qualitative (8 - 24x glycine is 3.2 - 9.6 mol/L, not a solution) | not run |

### 2.7 Fitted time courses

| source | what it shows | for / against | engine |
|---|---|---|---|
| Schieberle 2000 Table IV (DFA; 100 C, 0.5 M phosphate pH 5) | FFT 2.0 / 3.1 / 110 / 132 ug at 30 / 60 / 360 / 720 min; R0 fits FFT almost exactly, C1 is ~2 dex low at 30 min | signature (d): any slow supply must keep the FFT shape | **seen** (M1 D1) |
| B60's FFT rows (Hofmann 1998, Zhai 2023, Kang 2026 and others in `thiol_submodel.yml`) | the FFT levels B60 was fitted on | C1 costs 36 nats on B60's 52 rows | seen |

### 2.8 Earlier engine records on the same parts

| record | part | result |
|---|---|---|
| supply screen (prereg v2 sec. 0B; plug-in, no fits) | A, B, R-K, R-sink | S1, S2 (power law and Kallen-shaped) and S3 fail (b) and (c); the diagnostic S3 + exit passes (a) and (b) (loss ratio up to 12.3); no structure reaches (c): TTCA 18 mM with free reactants near 1 mM needs K_app >= ~19 L/mmol, 300x Cao's |
| M1, EX-061 | A (C1 vs R0) | "not separated", artefact-driven; evaluated exactly the reservoir favours C1; C1 3 - 5x too slow, loses the FFT shape |
| B58, B59b | A | Cao's K by detailed balance; de Roos's route; on Zhai 2019 the CAMA formation must be 2 - 5x slower (B59b prereg) |
| CF-076 (EX-050, EX-054, cluster 5, resolve2) | A, R-K | no single pH-flat exchange holds dilute and concentrated pots; Kallen-shaped K(pH) admitted provisionally (P 0.25), conditional on a low release Ea; catalysed exchange and a quench operator rejected |
| CF-063 | P1 - P4, R-cal | the deficit is 1 : 1 and formed within 5 min; P(adduct rather than a destructive sink) ~0.65 |
| CF-107 | all | whole-model refusals against instant supply are not evidence against a reservoir |

**Sources left out, with reasons.** The oxidant, O2 and disulfide parts (M1's C3, B11, B27) are not supply. Lipid,
buffer and hexose cells are M4 / M5 (the Cys + glucose barrier of Zhang 2026 enters only as `r_cama_form`'s existing
barrier). Storage and low-water cells (Kang's storage block, Zhai 2020, Zhai 2023c's dehydration, Mizutani's TEST
cell) are context. Dossiers outside F04, CF-063 / CF-076 and the files named above were not screened one by one: a
stated gap.

## 3. The rivals

| id | explanation | (b): TTCA loss vs label pH | (c): Zhang's missing ~19 mM | (a): Cao's K under drift | (d): FFT shape |
|---|---|---|---|---|---|
| **R0** | instant supply, the incumbent (B60) | no TTCA in Cys + Xyl pots; fed TTCA's loss pH-flat | no TTCA pool at all: reactants consumed in seconds | no K | fits (the reference) |
| **C1** | M1's reservoir, release free | pH-flat (ratio <= 1.01 in the screen) | TTCA lost when free Cys is low; no pool | flat by construction | lost at early times (2 dex at 30 min) unless the faster release restores it |
| **X** (claim) | C1 + the base-favoured exit (X-cama or X-deg) + P1 - P4 inside their bands | rises with pH; must also give the acceleration under Kang's acid drift | whatever P1 - P4 reach inside their bands (sec. 5.2 arithmetic: P1 - P3 at most ~7 mM cysteine and ~5 mM xylose, against 19) | flat if the exit is slow at pH 4.2 - 4.7 | open |
| **X0** | C1 + the exit, no pool | as X | not reached; Zhang's blank declared not representable | as X | as X |
| **R-cal** | X0 + an analytical cause: Zhang's free Cys and Xyl under-read at low level (CF-063's identical curve; the 0.223 mmol/L zero-signal point; ELSD non-linearity) | as X | no pool: at K = 0.060 L/mmol and TTCA 18.2 mM the free reactants should be ~17 mM each, which closes the balance to ~93 % (arithmetic) | as X | as X |
| **R-sink** | X0 + material destroyed, not held: S3's 1 : 1 drain on free Cys and PENT, and a cysteine-only sink (CF-076 explanation 5; Wei 2022) | as X | destroyed, so never returned on dilution | flat if the drain is slow at 90 C | open |
| **R-K** | Kallen-shaped K(pH) with the exit (CF-076 explanation 1 / 2) | K rising with pH stabilises TTCA; needs the exit to outrun it | a K ~300x Cao's at pH 6.5 holds TTCA at 18 mM against ~1 mM free reactants | moves <= 0.075 dex over Cao's drift for pKa_T >= 5.2 (screen arithmetic) | open |

X, R-cal and R-sink can hold together in parts; the case asks which parts the decisive data need. The stop rule is
the plan's: if no structure passes (a), (b) and (d) within ~25 h, supply is "not identifiable on disk". (c) is not in
the stop rule; it decides between X, R-cal and R-sink.

## 4. The decisive data

### 4.1 Why the screen's pots cannot decide

Cao 2017, Zhai 2021 and Zhang 2025b are **seen**: the screen evaluated the engine on all three across every
structure's band, and M1 scored Zhang's blank and Cao. Two of them also feed constants: Cao fixes K and anchors
`r_cama_form`, and Cao's drift and pool closure set the screen's m and drain ceilings. A comparison on them can only
say whether a structure is feasible (development grade), not which rival is right. The same holds for every other
supply measurement already run under the engine: Zhai 2019 Fig. 3A / 3B, Zhai 2023's Fig. 3a and Fig. 5, Kang's Fig.
S3 / S4, Schieberle 2000 (D1), B60's FFT rows, and de Roos 2005 Fig. 4 (spent as a constant).

### 4.2 Unseen data on disk: none is decisive

The held sources not yet run under the engine were screened for a measurement on which the rivals differ:

| source (not run) | what it could test | why it is not decisive |
|---|---|---|
| Hou 2017 (fed TTCA vs fed CAMA, 140 C, pH 5.5) | whether TTCA's thiols pass through CAMA (a lag in fed TTCA's MFT / FFT against fed CAMA's, same vial, same detector: a scale-cancelling ratio) | one pH, so silent on part B and on the pool; peak areas at 30 min only; the dossier holds no numbers (it needs keying from the PDF). A side unit for X-cama vs X-deg at most |
| Wei 2022 (GX-ARP + Cys, 110 - 130 C) | a cysteine-only sink (R-sink) against a 1 : 1 adduct | its sugar supply runs through the glycine Amadori compound's reversal, which the engine does not carry; pH not measured; unkeyed |
| Hou 2026 (TTCA and cystine after a 120 C step) | the cystine member's size | it sets P3's band (sec. 5.2), so it cannot decide |
| Kang 2026 end pH; Table S5; Table S6 | the pH trajectory; FFT / MFT vs pH; storage loss vs pH | the end pH is an input; S5 is compromised; S6 is storage |
| Zhai 2020; Zhai 2023c | TTCA -> CAMA vs pH (spray-drying); fed TTCA at pH 7.4 | outside the regime; no time course |

**No unseen decisive data exist on disk.** The plan gave P ~0.6 to that; after this screen of the held sources, P
~0.85 (the residual is the dossiers not screened one by one). **The verdict of M1c is therefore development-grade**:
a structure can be found feasible or infeasible on (a) - (d), and the rivals can be ranked on the development data,
but nothing on disk can ship a supply structure as decided. What would decide it is the designed experiment in 4.4.

### 4.3 The new papers (from the queue only, declared before reading)

| paper (queue description) | role, fixed now | the measurement that would make it decisive | which rivals it separates | P it holds that measurement |
|---|---|---|---|---|
| **belikov1967** (stability of L-cysteine + carbohydrate condensation products; possibly K for several sugar thiazolidines; title only) | **decisive candidate** | an equilibrium constant (or ring : open-chain ratio) for a cysteine-PENTOSE thiazolidine at >= 2 pH values, or >= 2 temperatures, in one study | R-K (K at pH 6.5 >= ~300x K at pH 4.4) against C1 / X / R-cal (K flat); a K(T) series also tests whether P1 grows with temperature | ~0.2 |
| **wlodek1993** (non-enzymatic hydrolysis of sugar TCAs; Cys-glucose among the least stable; released Cys + aldehyde -> a bis-Cys thioaminal pool) | **decisive candidate** | the release (hydrolysis) rate of a pentose TCA at >= 2 pH values between 6 and 8 at one temperature, cysteine measured | whether the base-favoured loss can be the release itself (Karpinski says no at 25 C). A bis-cysteine pool is 2 Cys : 1 aldehyde, so it cannot be Zhang's 1 : 1 deficit; that part is evidence only | ~0.35 |
| **davidek2003** (HPAEC-ECD / DAD of sugar + amino acid + Amadori together; detection floors) | **decisive candidate only if** it reports a cysteine + pentose pot; otherwise evidence for floors of its own assay | free cysteine, free pentose, TTCA / CAMA quantified together by an assay independent of ELSD, at >= 2 times | R-cal against X / R-sink: an independent assay that closes the balance favours R-cal; one that reproduces a 1 : 1 deficit favours X or R-sink | ~0.1 |
| gajda1990 (protonation constants of 2-(polyhydroxyalkyl)TCAs, pH 1.5 - 8, potentiometric) | **band source**: TTCA's ring-N pKa for R-K's Kallen shape (replaces Radomski's inferred 5.2 - 6.1) | -- | -- | -- |
| canle1996 (ring-opening k and ring / open K vs pH for N-aryl thiazolidines, buffer catalysis, 25 C) | **band source**: the class shape of P1's ring-chain K and of ring opening | -- | -- | -- |
| szilagyi1979 (C-2 inversion and C-4 epimerisation in TCAs; ring-chain rates by NMR) | **band source**: P1's exchange speed; P2's epimer share | -- | -- | -- |
| davis1991, davis1991b (penicillin thiazolidine ring opening: pH-rate profiles) | **band source**: the pH shape of ring opening (acid-catalysed below ~5, flat 6 - 12 per the queue) | -- | -- | -- |
| carta1996 (L-cystine solubility vs pH and NaCl) | **band source**: P3's ceiling and whether cystine precipitates on cooling (an observation-operator question) | -- | -- | -- |
| sohn1995 (ammonia from amino acids heated in water) | **band source**: R-sink's cysteine-only arm (cysteine deamination); also CF-084 (M1b) | -- | -- | -- |
| davidek2004 (pentose-Amadori assay; xylose + glycine, 90 C, pH 6, ARP ~12 mol % at 15 min) | **evidence only** (glycine, not cysteine): neither decisive nor a band | -- | -- | -- |

P(at least one of the three candidates supplies a decisive measurement) ~0.5. None of the eleven sets a band in
this case: a band source enters only through an amendment registered before any engine result of M1c exists.

### 4.4 The designed experiment (W6)

One laboratory, three arms, n = 3 independent heatings per cell, sealed vessels with stated fill and headspace (N2 and
air arms in Arm 2), pH measured on the cooled scale for every sample and in situ where an electrode allows, every
assay's LOD / LOQ reported.

| arm | conditions | species | time points | decides |
|---|---|---|---|---|
| **1. the exit's pH law at held pH** | purified TTCA (equilibrium epimer mixture, NMR-characterised) 10 mmol/L; pH 4.5, 5.5, 6.5, 7.5, 8.0 held by pH-stat, or in phosphate at 0.05 and 0.2 M to extrapolate the buffer term to zero; 100 and 120 C (90 C at pH 4.5 to bridge Cao); +- 10 mmol/L [13C5]xylose | TTCA (each epimer), CAMA (each isomer), free Cys derivatised at the quench (N-ethylmaleimide, isotope-labelled internal standard) so it cannot re-form TTCA, free xylose (HPAEC-PAD), [13C5]-TTCA (release + recombination), cystine, H2S (trap), 1-DX / 3-DX (OPD on the quenched sample), A420, pH | 0, 2, 5, 10, 20, 40, 60, 90, 120 min | part B at fixed pH (no drift): its slope n; which step carries it (CAMA appears: X-cama; deoxyosones or browning without CAMA: X-deg; free Cys appears and stays: the release, which the analogues say is not base-favoured); the acceleration (an induction under constant pH means a pool such as P1) |
| **2. Zhang's pool in a second laboratory** | Cys + Xyl 40 + 40 mmol/L, 0.2 M phosphate pH 6.5, 120 C, 5 mL in 15 mL (Zhang's design); the same recipe at 0.4 M and 90 C (Cao's design) as a bridge; a [13C5]xylose run | free Cys (as in Arm 1), free xylose (HPAEC-PAD, an assay independent of ELSD), TTCA and CAMA by LC-MS/MS with labelled standards, total releasable cysteine (a validated exchange that opens thiazolidines, checked on TTCA standards), cystine in solution AND in any pellet, 13C NMR of [13C5]xylose samples quenched at 10 and 50 min (open-chain and epimer signals); a 20x dilution of the 50-min sample at 120 C followed for 30 min; immediate cold dilution vs 4 C storage 0 / 4 / 24 h | 0, 1, 2, 5, 10, 20, 30, 50 min | (c): R-cal (an independent assay closes the balance) vs X (a 1 : 1 pool that NMR sees and that gives Cys and Xyl back on dilution) vs R-sink (deficit not returned on dilution, not releasable) |
| **3. K(pH) and K(T)** | Cys + Xyl equilibrated from both sides (from TTCA, and from Cys + Xyl), pH 4.5, 5.5, 6.5, 7.5, at 25, 60 and 90 C; TTCA ring-N pKa titrated at 25 and 90 C (CF-076's cheapest datum) | all four species and pH at temperature | to a plateau (checked by two successive samples) | R-K vs a flat K; P1's temperature law |

Arms 1 and 2 decide the structure; Arm 3 sets K's pH law (a band source for the next case if Arms 1 - 2 are its
decisive data, or decisive for R-K if Arm 3 is the only arm run). Request it through W6, with CF-076's existing
designed experiment (TTCA 10 vs 100 mM, +-phosphate, +-morpholine) folded into Arm 1.

## 5. Tolerances from the data

### 5.1 The measurement model per unit

Rules (M1.1): tau (H0's tau_ABS 0.964 dex) only on absolute levels compared across laboratories or against an
engine level; omega 0.390 dex = 0.898 ln (EX-059) in every unit's primary covariance, per printed value; floors from
each assay's LOD / LOQ or calibration range, never a data point; exact or Monte-Carlo densities (PSIS k-hat < 0.7, ESS
>= 100) on every form that can flip a sign.

| unit | value basis | reported noise | error terms in the primary | floor | note |
|---|---|---|---|---|---|
| (a) Cao 2017, log K_app vs measured cooled pH, 0.5 - 5 h | quotient of three absolute levels; its slope on pH is level-free (a multiplicative calibration bias cancels in the slope) | n = 3, RSD < 10 % per species (0.095 ln); read error +-0.004 mol/L (Cys, Xyl), +-0.005 (TTCA), +-0.03 pH | omega on each of the three printed values: log10 K sd 0.675 dex per time point; no tau | none needed (Cys 0.056 - 0.076 inside its 0.01 - 0.1 mol/L calibration) | **under omega the slope's sd is ~1.5 dex / pH** (0.675 / sqrt(0.19), the pH deviations' sum of squares): (a) as a likelihood unit has almost no power; the m <= 0.18 bound holds only on measurement noise. Whether omega enters per printed value or per sample (shared across species) is open question 1 |
| (b) Zhai 2021 Fig. 4b, ln(C(t) / C0), pH 5.5 / 6 / 7 / 8, +- Xyl | C / C0 with C0 the plotted 10.00 mmol/L charge: calibration cancels; no tau | SD (n = 3) printed for TTCA alone at pH 5.5 / 6 / 7 only; for pH 8 and the +Xyl series a transfer, the median printed relative SD of the same figure (~3 - 4 %), sensitivity x0.5 / x2; read error +-0.005 mmol/L | omega per point | TTCA LOD 2.21 mg/L = 0.0087 mmol/L (Zhai's printed); "0.00" points are censored below it | arithmetic on the measured series alone: a pH-flat rival (the pH 5.5 series as its stand-in) misses the pH 7 and 8 points by ~6 nats under omega, so (b) keeps its power. The pH 7 points at 100 / 120 min print SD 0 (not measured as SD): the transfer applies |
| (c) Zhang 2025b: TTCA, CAMA, free Cys, free Xyl, 5 - 50 min | free Cys / Xyl: ln(C(t) / C(5 min)), first-point-normalised (prereg v2); TTCA and CAMA: absolute levels within the lab | n = 3, error bars small and not transcribed; read error +-0.3 mmol/L (TTCA), +-0.05 (CAMA), ln read error on 2g / 2h | omega per point; tau on TTCA and CAMA levels when compared with an engine level (absolute basis) | Cys / Xyl: Zhai 2021's LOD 0.025 / 0.020 mmol/L as cited by the prereg v2 (UNVERIFIED: the zhai2021 dossier prints only TTCA's), alternative 0.223 mmol/L (the printed ELSD curve's zero-signal point); TTCA / CAMA: the calibration floors 0.1 mmol/L | the K gap: Zhang's K_app 18.6 L/mmol against Cao's 0.060 is 2.49 dex; same laboratory, so omega, not tau: sd sqrt(2) x 0.675 = 0.96 dex, **~2.6 sd**: real, not overwhelming, and one pot |
| (d) Schieberle 2000 FFT and MFT ratios 360 / 30, 720 / 360 | ratios within one series: no tau; bivariate (shared 360-min point) | no SD or n printed; M1's 0.094 ln per point (read error) | omega per point | 1 ug/L (the lab's "< 0.1 ug" per 100 mL) | as M1 and the prereg v2 |
| (d) B60's FFT rows | absolute levels, several labs | each row's printed SD or its declared transfer | omega + tau (absolute) | each row's stated floor | the in-sample guard |
| decisive, Arm 1 | ln(C(t) / C0) per species and pH; the slope n from the within-study pH contrast | the lab's replicate SD per point | omega per point; no tau | the arm's own LOQs; censored likelihood below them | exact densities |
| decisive, Arm 2 | the closure fraction (sum of species / charge) at each time; ln(C after dilution / C before) | replicate SD; recovery of each assay from spiked standards (reported) | omega per printed value; a per-assay calibration term from the spike recoveries (the closure sums absolute levels from different calibrations); tau only if compared with Zhang's own levels | the arm's LOQs | isotope-dilution assays keep the calibration term small |
| decisive, Arm 3 / belikov1967 | ln K ratio between two pH values (or temperatures) inside one study | as printed; if none, a declared transfer (Cao's RSD < 10 % per species), sensitivity x0.5 / x2 | 2 omega^2 per ratio of two K (each K three values); no tau | the paper's floors | -- |
| decisive, wlodek1993 | ln ratio of release rates between pH values inside one study | as printed; else the transfer above | omega on the underlying levels; no tau | the paper's floors | 37 C: transferred by the exit's declared Ea band only as a sensitivity |
| decisive, davidek2003 (if a Cys pot) | closure fraction within one study | its printed RSDs and LODs | omega per value; per-assay term from its printed recoveries | its printed LODs | -- |

Across laboratories (any unit that compares a level from one lab with a constant from another): H0's tau_ABS.

### 5.2 The bands and their sources

No band below comes from Cao 2017, Zhai 2021 or Zhang 2025b (screened) or from any candidate named in 4.3
(decisive). Where the only source is a screened pot, the value enters as a fixed measured constant (and the case says
which signature it then cannot test) or as a required-agreement sensitivity, never as a band.

| coordinate | band | source | kind |
|---|---|---|---|
| n, the exit's pH slope (dex per pH unit, cooled scale) | **0.37 - 0.92** | Martins & van Boekel 2003 Part II Table 3: per-step slopes of amine-sugar Amadori steps, 100 - 120 C (`trunk_conditions.PH_SOURCE`) | same-class (proton-transfer steps of an amine-sugar adduct); sensitivity 0 - 1 (one hydroxide, mechanism) |
| the exit's level at pH 5.5 | X-cama: `r_cama_form`'s shipped constant (3.29e-3 /min at 90 C, Ea 100.9 kJ/mol), fixed; X-deg: `k_ttca_deg` B60 +- 1 dex | B59b anchor on Cao's CAMA rise (measured); Zhang 2026's Cys + glucose Amadori barrier (hexose transfer); X-deg's width has no source | X-cama: measured constant from a screened pot, so (a) tests only K's flatness, not the CAMA rise (declared); X-deg: numerical |
| `k_ttca_cys`, the release | log10 k(145 C) = B60 +- 3 dex; Cao's floor (>= 7e-3 /min at 90 C) as a required-agreement sensitivity only | stage-1 convention; Cao arithmetic 4 | numerical (primary), screened (sensitivity) |
| K | 0.060 L/mmol fixed, pH- and T-flat in C1 / X / X0 / R-cal / R-sink; 0.045 - 0.075 as sensitivity | Cao 2017, the only cysteine-pentose K on disk | measured constant from a screened pot: (a) then tests flatness under drift, not K's level |
| R-K: pKa_T in the Kallen shape | 5.2 - 6.1 | Radomski 1989 (room-temperature charge-balance inference, as CF-076); gajda1990 replaces it by amendment if read | inferred class |
| P1: f_oc = [TTCA_oc] / [TTCA] | 0 - 0.05 at <= 37 C; temperature law Delta H 0 - 60 kJ/mol | Radomski 1989 (acyclic form not seen by NMR, 37 C); Karpinski 1990 (no open-ring thiol, LOD ~1 uM, 25 C) | same-class at <= 37 C; **numerical** in temperature (no thiazolidine ring-chain T series exists) |
| P2: epimer share not integrated | 0 - 0.20 of TTCA (2R : 2S = 4.06) | Radomski 1989, 37 C | measured (TTCA itself); operator switch on / off |
| P3: Cys -> 1/2 cystine | first-order ceiling 7e-4 /min at 120 C (<= 8 % of Cys in 120 min); barrier 40 - 135 kJ/mol | Hou 2026 L-MRPs (Jiangnan, 120 C step); `sugar_trunk` class for the barrier | same-lab-network ceiling; numerical barrier |
| P4: second CAMA isomer; `r_tdp_cys_trap` | 0 - 1 of CAMA; `k_tdp_cys_trap` 0 - its EX-046 value | none measured | **numerical**: a verdict that needs P4 does not ship |
| R-cal: under-read factor on Zhang's free Cys / Xyl below 5 mmol/L | 1 - 30 | none on disk (hu2012, queued and unread, is a candidate class source by amendment) | **numerical** |
| R-sink: 1 : 1 drain on free Cys and PENT | 1e-6 - 3e-3 /min at 90 C; barrier 40 - 135 kJ/mol | numerical; Cao's closure ceiling as a required-agreement check | numerical |
| R-sink: cysteine-only sink | 0 - 4e-4 /min at 90 C | Cao's closure (screened): sensitivity only; sohn1995 by amendment | numerical (primary) |
| pH trajectory, Zhai-design pots | linear in time from the label pH to Kang 2026's measured end pH at 120 min (5.1 / 4.9 / 4.7 for 5.5 / 7 / 8; 5.0 interpolated for 6), cooled scale; sensitivity: the engine's charge-balance pH | Kang 2026 p. 3242 (same laboratory and design; its temperature unresolved, 100 or 120 C) | measured end points, declared shape |
| pH trajectory, Cao 2017 | the keyed cooled readings (Fig. 4c) | Cao 2017 | measured |
| pH, Zhang 2025b | label 6.5 (0.2 M phosphate; no reading printed); sensitivity: the engine's charge-balance pH | -- | declared |

**What the bands allow (arithmetic, no engine).** Inside the non-numerical bands, P1 holds at most ~0.9 mM (5 % of
18 mM TTCA at <= 37 C), P2 at most ~4.5 mM, P3 at most ~1.4 mM of cysteine by 50 min (7e-4 /min) and no xylose: ~7 mM
of cysteine and ~5 mM of xylose against the 19 mM of each that (c) needs. So the claim reaches (c) only through a numerical band (P1 at
high temperature, P4), or not at all. P(X passes (c) inside non-numerical bands) ~0.1; R-cal's single bias closes the
balance by construction, but its band is numerical too. Either way (c) cannot ship on disk: Arm 2 decides it.

## Open questions for the build

1. **Omega per printed value or per sample?** Under M1.1 as written, signature (a) has almost no power (sec. 5.1).
   If model error is shared across species within one sample, K's flatness keeps power. Declare before the screen.
2. **The cooled scale.** n is banded on cooled-scale pH (Martins, Cao and Kang report that scale), but the chemistry
   acts at the in-situ pH. The primary uses the engine's cooled-scale pH of the current composition
   (`observations.ph_cooled`); the in-situ argument is a required-agreement sensitivity. Is that acceptable?
3. **X-cama against de Roos and Zhai 2020.** De Roos says the Amadori route is favoured below pH 7, Zhai 2020 says
   base favours TTCA -> CAMA (dehydration). Both variants run; neither source chooses.
4. **The acceleration.** Zhai 2021's curves speed up while Kang's pots drift acid. If no variant of part B produces
   the acceleration on the declared trajectory, is (b) met by the 80-min ratio alone (the screen's form) or must the
   shape pass too? Default: the 80-min ratio decides feasibility; the shape is reported.
5. **Cao as anchor and signature.** K and `r_cama_form` come from Cao's blank, which is also signature (a). Sec. 5.2
   says what (a) then cannot test. A reviewer may prefer to drop (a) as a test; the plan's stop rule then rests on
   (b) and (d).
6. **The FFT shape (d) and the release.** C1 lost the FFT shape at 30 min; the M1 review attributes it to a release
   3 - 5x too slow. If freeing `k_ttca_cys` restores (d) in C1 alone, part B is still needed for (b), but the case's
   claim then rests on (b) only.
7. **Keying Hou 2017.** Its fed TTCA vs fed CAMA ratio is the only rank-2 datum on disk not yet run. Key it (2 - 3 h,
   two readers) as a side unit for X-cama vs X-deg before the screen, or leave it out with this reason.
