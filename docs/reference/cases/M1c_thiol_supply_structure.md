# Case M1c: the thiol supply structure (a base-favoured TTCA exit and Zhang 2025b's missing pool)

**FROZEN 2026-11-05 (items 1–5), before any of the 11 new papers was read.**

*Plan step M1c ([STATUS_AND_PLAN sec. 5, W2](../../guides/STATUS_AND_PLAN.md), as re-sequenced 2026-11-05; error E10,
with E12 and CF-076 / CF-063 / CF-107). Items 1 - 5 of sec. 7, plus the bands every new coordinate carries, written
2026-11-05 before any result of the change exists and before any newly downloaded thiol-supply paper was read. No
engine was run for this case, on any pot. Items 6 - 9 are appended later and never edited into this part.*

**Version 2 (2026-11-05), after an independent blind review (FREEZE AFTER FIXES) of v1 (88457d7e).** Still before any
result and under the same blinding. v2 adds Zhao 2017b and nine other held sources, and a screen of all 89 held
dossiers that mention thiazolidine or TTCA (sec. 2.9). It re-specifies R-cal and adds R-rel. It re-derives the band
on n and re-references the X-cama anchor. It changes the primary forms of (b) and (c), demotes (a) to a consistency
check, and restates the stop rule. It extends the designed experiment and moves omega to 0.419 dex. The revision log
is at the end.

**Blinding.** Eleven newly downloaded PDFs bear on thiol supply: `canle1996`, `gajda1990`, `belikov1967`,
`wlodek1993`, `davidek2003`, `davidek2004`, `carta1996`, `szilagyi1979`, `davis1991`, `davis1991b` and `sohn1995`.
None was opened, and no extraction of them was read. What this case says about them comes from the private download
queue only. The queue gives titles for all of them, abstract-level descriptions for some, and the title only for
`belikov1967`, `carta1996`, `szilagyi1979` and `sohn1995`. Sec. 4.3 fixes now which of them could be decisive, and on
what measurement, and which are band sources. A paper named as decisive can never set a band, and a band source can
never decide.

**Seen by the author.**
- **Engine results:**
  - the supply screen's outcomes as the parked M1 + M2 prereg v2 records them (sec. 0B table: each structure's
    feasibility on signatures (a) - (c));
  - engine predictions on Zhang 2025b's blank pot (best TTCA 11.6 / 8.8 mmol/L under S1; free Cys >= 3.2 / ~10 / 15
    mmol/L when TTCA and CAMA pass under S2 / Kallen / S3);
  - the engine's cooled-scale pH in Cao 2017's pot (~0.2 low at 30 min under C1);
  - M1's decisive comparison and results review (EX-061: R0 and C1 per-species residuals on Zhang 2025b and
    Schieberle 2000);
  - the no-refit probe on Zhai 2019 Fig. 3A (TTCA -0.29 dex mean, CAMA +0.71 dex) and B59b's prereg check (CAMA
    formation 2 - 5x too fast), relayed from the repository's records by a reading agent;
  - for Zhao 2017b, C1's predictions (TTCA 252 vs 355 mmol/L at 30 min, 120 vs 290 at 300 min), relayed by a reading
    agent from M1's fit and side-effect records.
- **Measured values only (no engine number):** the dossiers named in sec. 2, the reading agents' screen lines for the
  89 dossiers of sec. 2.9, the F04 atlas file, and the queue's 2026-11-05 section.
- **Not seen:** any engine result on Hou 2017, Wei 2022, Hou 2026, Zhai 2020, Kang 2026's storage and pH-ladder
  tables, or the eight analogue and sink dossiers added in v2. Also not seen: Mizutani 1960's F04#26, E10's only v2
  TEST contrast, which is read in item 7 only. The screen's own record (`kinetic_core_supply_structure_signature_screen.md`)
  did not exist when this was written.

## 1. The claim

In a heated cysteine + pentose pot, cysteine is held in a reservoir of TTCA and the cysteine Amadori compound (CAMA),
as in M1's C1. The claim adds two things.
1. **A base-favoured exit.** TTCA leaves the reservoir through an irreversible exit that is faster at higher pH, so
   fed TTCA is lost ~10x faster at label pH 8 than at 5.5 (Zhai 2021 Fig. 4b).
2. **An account of Zhang's missing pool.** Zhang 2025b's ~19 mmol/L of each reactant that no measured pool holds is
   accounted for by named pool members. Each member enters with a measured or same-class constant or a declared
   band, never as one free pool sized to the gap. The size of the account is an output. If the members' bands cannot
   reach the gap, the claim says so rather than widening them.

Network steps that change (names from [sulfur.py](../../../src/kinetic_core/sulfur.py) and
[parameters_sulfur.py](../../../src/kinetic_core/parameters_sulfur.py)):

| part | step(s) | shipped (R0, B60) | C1 (M1) | M1c |
|---|---|---|---|---|
| A. the reservoir | `r_ttca_form` / `r_ttca_cys` (Cao's K, detailed balance); `r_cama_form`; `r_cama_dpo` / `r_cama_tdp` | off (`TTCA_FORMATION`, `CYS_AMADORI_ROUTE` disabled): instant supply through `r_pent_dpo` / `r_pent_tdp` | on; K 0.060 L/mmol, pH- and T-flat; `r_cama_form` pH-flat (3.29e-3 /min at 90 C, anchored on Cao's blank, Ea 100.9 kJ/mol); CAMA decay on de Roos 2005 Fig. 4 | as C1, with the release `k_ttca_cys` free and K(T) banded (sec. 5.2) |
| B. the base-favoured exit | a factor 10^(n (pH_c - pH_ref)) on ONE irreversible TTCA exit. **X-cama**: on `r_cama_form`, pH_ref 4.4 (Cao's time-weighted cooled pH, where its constant was anchored). **X-deg**: on `r_ttca_deg`, pH_ref 5.5 | none | none | both variants built, both run, no choice between them on scores; pH_c = the engine's cooled-scale pH (sec. 5.2) |
| C. reversible members (1 : 1 Cys : Xyl) | **P1a** open-chain / Schiff-base TTCA; **P1b** the cysteine N-xylosylamine (a TTCA isomer, MW 253; it would give the same 254 > 122 transition, so it is invisible only if it elutes apart from TTCA); **P2** the minor C-2 epimer (observation operator only); **P4a** a second CAMA isomer outside Zhang's MRM transition | absent | absent | each inside its band; they return Cys and Xyl on dilution |
| D. irreversible members | **P3** Cys -> 1/2 cystine (cysteine only); **P4b** 3-deoxypentosone + Cys (`r_tdp_cys_trap`, 1 Cys : 1 C5) | absent / zero | absent / zero | each inside its band; they return nothing on dilution |
| not changed | Cao's K level at 90 C; the CAMA decay law; the thiol loss law (`k_mft_decay` as M1 left it); the pH balance (M1b(i)) | | | |

The run uses the measured pH trajectory where one exists (Cao 2017; for the Zhai-design pots, Kang 2026's end points
with three declared shapes, sec. 5.2), whatever M1b(i) finds. Every pH comparison states its scale, cooled or in situ.

## 2. All the evidence, by kind

"Engine" says whether the measurement has been run under the engine (**seen**) or not (**not run**). Laboratories are
as in [lab_ids.yml](../contrast_benchmark/lab_ids.yml). For E10 the frozen v2 split
([split_v2_per_change.yml](../contrast_benchmark/split_v2_per_change.yml)) has three DEV labs: BTBU (Xie), the Ho
network (Jiangnan) and Kallen (Pennsylvania). It has one TEST contrast, F04#26 (Mizutani 1960, Hokkaido), which is
used in item 7 only.

### 2.1 Isotope labelling

| source | what it shows | for / against | engine |
|---|---|---|---|
| Zhai 2021 Fig. 7a (13C5-xylose + fed TTCA; CF-091) | label crosses into TTCA, so release and recombination occur at 100 C. Quantified only as a difference of two pots | for a reversible release (C1, X); at face value it disfavours CF-076's E4c point (z 6.3) | seen (CF-076 cluster 5) |
| Zhai 2023 Food Chem Table 2 (13C5-xylose + fed TTCA, 120 min) | exogenous xylose carbon in thiophenes: 0 % (100 C), 19 - 21 % (120 C), 34 - 45 % (140 C); no n | release + recapture grows with T | not run |
| Zhao 2019 (BTBU; [13C5]-TTCA, 140 C, pH 5.5, 30 min) | 97 % of MFT carries the intact TTCA C5 | TTCA's C5 reaches MFT (C1, X) | not run |
| none | no label on the base-favoured loss's products or on the missing pool | gap | -- |

### 2.2 Fed intermediates

| source | what it shows | for / against | engine |
|---|---|---|---|
| **Zhai 2021 Fig. 4b** (Jiangnan; TTCA 10 mM +- Xyl, unbuffered, label pH 5.5 / 6 / 7 / 8, 100 C) | fraction of TTCA left at 80 min: 0.819 (pH 5.5) vs 0.147 (pH 8). First-order equivalent over 0 - 40 min: 0.0017 -> 0.0164 /min (~0.4 dex / pH). Every curve ACCELERATES (pH 7: 0.0045 /min over 0 - 20 min, 0.049 over 80 - 100). SD printed only for TTCA alone at pH 5.5 / 6 / 7 | for part B. The acceleration is **against a first-order exit on a pH that drifts acid** (that would decelerate). It is for an induction: a P1-type pool that builds up, or autocatalysis | **seen** (signature (b); B16 / route rows; CF-076) |
| Zhai 2021 Fig. 5 (same pots; OPD, 3-DX / 1-DX on a diacetyl surrogate, semi-quantitative) | standing 3-DX per TTCA lost at 120 min: 19 / 17 / 5 / 0.8 % (pH 5.5 / 6 / 7 / 8). 1-DX at most 106 umol/L. A420 doubles from pH 5.5 to 8 | the fast base-favoured loss goes to neither deoxyosone in measurable amounts: into browning or an unmeasured product | seen (deoxyosone rows); a declared routing check (sec. 5.1) |
| Zhai 2021 Fig. 3a (pH 5.5; 100 / 120 / 140 C) | zero-order k 0.0274 / 0.0651 / 0.0813 mmol/L/min. The printed Ea of 80.99 kJ/mol is not reproduced (~35 kJ/mol) | the exit's temperature law at one pH | seen |
| Kang 2026 p. 3242 (Jiangnan; TTCA 10 mM, the Zhai design) | **measured** final pH 5.1 / 4.9 / 4.7 from label 5.5 / 7 / 8, cooled. The series temperature is unresolved (100 or 120 C) | the label-pH contrast is a contrast in the early minutes only | not run (an input, sec. 5.2) |
| Kang 2026 SI Fig. S3 / S4 | free Cys released from fed TTCA peaks at <= 16.3 mol % (100 - 140 C). The F04 atlas records S3 as the same data as Zhai 2021 Fig. 6 | for a release that recombines | seen (B2.1 fits; `k_ttca_deg`'s B60 value rests on it) |
| Kang 2026 SI Table S5 | FFT / MFT 2.91 / 3.02 / 2.73 over label pH 5.5 - 8 | weak: S5's pH-7 column is S4's 100 C column | not run; compromised |
| Hou 2017 (BTBU; fed TTCA, CAMA, Gly-ARP, 30 mM each, 0.2 M phosphate pH 5.5, 140 C, 5 - 30 min) | reactivity Gly-ARP > CAMA > TTCA. TTCA alone and CAMA alone both make MFT and FFT (peak areas) | for TTCA's thiols passing through a slower step than CAMA's; one pH | not run; unkeyed |
| Zhang 2023b (BTBU; fed TTCA + Met, 0.01 M phosphate, initial pH 3 - 8, 115 C) | volatiles only; no TTCA, Cys or H2S measured; levels disagree between figures by up to ~700x | none usable for supply | not run |
| de Roos 2005 Fig. 4 (Givaudan; fed CAMA, boiling 0.5 M phosphate, pH 3 - 6) | CAMA decay U-shaped: slowest at pH 5 (4.5e-3 /min), 1.43e-2 /min at pH 6 | already the engine's CAMA decay law | spent (B59b constant; keyed held-out for E10) |
| de Roos 2005 text | "At pH < 7, the route to Amadori compounds is preferred"; the thiazolidine is favoured above pH 9 | **against X-cama**, qualitatively | not run |
| Zhai 2020 (Jiangnan; spray-dried TTCA intermediate) | TTCA -> CAMA conversion favoured at pH 7.5 - 9.5 (CAMA 6 % -> 47 % at pH 9.5) | **for X-cama**, but at low water during dehydration | not run |
| Zhai 2019 Fig. 3B; Zhai 2023c | dicarbonyls only; no TTCA time course | context | seen; not run |
| Wei 2022 (Jiangnan; GX-ARP 5 + Cys 5 mM, water, label pH 7, 110 - 130 C) | at 120 C and 100 min, 4.24 mM Cys consumed and 0.48 mM TTCA formed (11 %). Free Cys 5 -> ~1.4 mM in 20 min | most consumed cysteine sits in no measured pool where little free xylose exists: for a cysteine-only sink (R-sink's cysteine arm) | not run; unkeyed |
| Zhou 2022 (Jiangnan; Cys 20 mM alone, 120 C, label pH 7, unbuffered) | 94 / 86 / 79 / 74 % of the Cys left at 30 / 60 / 90 / 120 min | the cysteine-only sink in water: ~2.5e-3 /min at 120 C | not run (only Zhou 2022's ARP + Cys rows are in the engine records) |
| Liu 2022 (Cys 40 + MGO 20 mM, 120 C, label pH 7.5) | 26.7 mM Cys lost between 15 and 30 min against ~1.4 mM MGO | a fast cysteine + dicarbonyl sink | not run |

### 2.3 Trapping

| source | what it shows | for / against | engine |
|---|---|---|---|
| Schubert 1936 / 1939 (iodoacetate, I2, room temperature) | an iodoacetate trap frees cysteine from the benzaldehyde thiazolidine (benzaldehyde separates after ~0.5 h; no ring left overnight). The xylose and arabinose rings take up I2 as fast as free Cys | the open chain is reachable in minutes at room temperature, so **a thiol-trapping quench can drain TTCA** (Arm 1's NEM control). No rate or equilibrium constant is printed | not run; no usable constant |
| Zhai 2019 Flag 1 (120 C, TTCA + Cys) | "TTCA disappeared completely" within 10 min, and a xylose peak appeared | against plain mass action; one reading is thiol addition to the open-chain Schiff base (P1a) | not run; qualitative |
| none | no in-situ trap of free cysteine, open-chain forms or cystine at process temperature | gap | -- |

### 2.4 Measured intermediates and lifetimes

| source | what it shows | for / against | engine |
|---|---|---|---|
| **Cao 2017** (BTBU; Cys + Xyl 0.4 M each, 0.2 M phosphate, 90 C, 0.5 - 5 h; TTCA by LC-MS SIM) | K_app 58.6 / 64.6 / 60.2 / 59.0 / 58.5 / 61.1 L/mol while the cooled pH drifts 4.67 -> 4.15. Slope of log K on pH +0.013 dex / pH, 90 % CI -0.155 to +0.180, on measurement noise alone. The pool over-closes by 6 - 17 % | a pre-equilibrium reservoir; K flat under drift | **seen** (signature (a); B59b anchor; M1) |
| **Zhao 2017b** (BTBU; the same recipe, vial and temperature; TTCA and CAMA by HPLC-ELSD; pH of a 6-fold dilution at 25 C) | K_app 59.2 / 60.3 / 62.6 / 66.8 / 71.2 / 78.0 L/mol while that pH falls 4.51 -> 4.10. Recomputed slope of log10 K on pH: **-0.290 dex / pH (90 % CI -0.331 to -0.248, n = 6, t4)**, outside Cao's CI. Pool +9 to +16 %; cysteine-loss ceiling <= ~2e-4 /min; TTCA -> CAMA net floor ~7e-4 /min; CAMA / TTCA 0.32 at 5 h (Cao 0.56) | **against "flat K"**: the same pot, run again with a different assay, gives a K that rises as pH falls. The opposite direction to R-K's Kallen shape. The dossier puts P 0.7 on separate physical runs | **seen** (M1 fit rows and side effects) |
| **Zhang 2025b** (BTBU; Cys + Xyl 40 + 40 mM, 0.2 M phosphate pH 6.5, 120 C, 0 - 50 min) | see the bullets below this table | for "a 1 : 1 pool or an analytical cause"; against cystine or H2S alone (cysteine only) and against sugar loss alone (xylose only) | **seen** (signature (c); M1; the screen) |
| Zhai 2019 Fig. 3A (Jiangnan; Xyl 165 + Cys 82.6 mM, water, label pH 7.4, 70 - 130 C) | TTCA 14.5 mM before heating. TTCA drains above ~110 C (130 C: 45.7 -> 6.9 mM). NMR finds no open-chain CAMA | the reservoir's T-law | seen |
| Zhai 2023 Food Chem Fig. 3a; Zhai 2023 JAFC Fig. 5 | regenerative free Cys from fed TTCA plateaus at ~1.0 mM (~10 %). Fed-TTCA CAMA peaks at 1.93 mM (120 C) | release small, CAMA small in water | seen |
| Radomski 1989 (Warsaw; TTCA among 15 sugar rings, water, 37 C) | 2R : 2S = 4.06. C-2 epimerisation relaxes in minutes (k ~0.20 /min, glyceraldehyde ring). The acyclic form is not seen by NMR | P1a small at <= 37 C; P2 <= 20 % of TTCA | band source |
| Roberts 1987 (ribose / xylose-Cys, ~25 C) | mutarotation fast enough to force a 3-min polarimetry convention; pyridine shifts toward the ring | qualitative; **no usable constant** | not run |
| Nagasawa 1981 (penicillamine + acetaldehyde, room temperature) | C-2 epimerisation takes hours (equilibrium within ~24 h) in neutral D2O. It takes < 1 h in 1 N HCl, so it is acid-catalysed, and it stops in 2 % Na2CO3. No free amino thiol on TLC | **against a base-favoured release**; a class bound for the P1 / P2 exchange (>= ~2e-3 /min neutral) | not run |
| Karpinski 1990 (Warsaw; voltammetry, deaerated, 25 C) | L-arabinose ring (AThz) release: 1.8e-4 /min at pH 6.0, 4.8e-5 at pH 7.4 (RSD up to 30 %). No open-ring thiol detected (LOD ~1 uM) | **release is not base-favoured above pH 6** | band anchor in CF-076's resolve2 |
| Bustillo Trueba 2018 (Leuven; Cys + aryl and Strecker aldehyde rings, 20 C) | ring loss bell-shaped in pH, maximum 4.4 - 5.2, nearly stable at 9.0 | as Karpinski | not run |
| Pesek 1975 (aliphatic TCAs, D2O, T not stated) | stable for days in neutral D2O and in 8 M HCl; decompose only in > 3 M NaOH | as above | not run |
| Fife 1991 (2-aryl thiazolidines, 30 - 90 C) | ring opening is fast and reversible; hydrolysis of the iminium is rate-determining and general-base catalysed, phosphate included | the form of a release law | not run |
| Kallen 1971 / 1971b / 1973 (Pennsylvania; Cys + formaldehyde, 25 C) | K_app rises ~20x from pH 5.6 to 7.5. k_r ~1e-6 /s. Phosphate catalysis saturates at x3.5 - 5 | the shape R-K borrows | seen (CF-076 battery) |
| Kim 1986 (Cys + cinnamaldehyde step 2, pH 8.0, 13 - 38 C) | K2 1.49e6 /M; Delta H2 +1.18 kcal/mol (+4.9 kJ/mol); one pH | a class anchor for K(T) only (sec. 5.2) | not run |
| Lienhard 1966 (thiol + aldehyde hemithioacetals, 25 C) | K_hemi 32 - 60 /M (acetaldehyde), 15 - 18 (isobutyraldehyde); a fast pre-equilibrium | **arithmetic excludes it as Zhang's pool**: at ~1 mM free Cys and Xyl, a 1 : 1 hemithioacetal holds ~K x 1e-6 M, i.e. < 0.1 uM, not 19 mM, even with K x 100 for a sugar carbonyl | not run |
| Wedzicha 1991 (3-deoxyhexosulose + Cys, acetate pH 5.5, 10 - 55 C) | binding is 1 : 1 (DH : Cys 1 : 0.96) and fast: initial rate 3.8e-4 M/min at 10 C, 24 mM DH. About 9 mM Cys is never liberated by 25 h | **a same-class band for P4b**: binding is not rate-limiting, so P4b's size is set by deoxyosone supply | not run |
| Edwards 1992; Wondrak 2002; Pripis-Nicolau 2000 | Cys + 3-deoxyosone adduct 74 % in 0.5 h (55 C, pH 5.5); Cys + phenylglyoxal 0.63 /M/s (37 C, pH 7.4); Cys + MGO ~87 % gone in 10 min (25 C, pH 3.5) | the same class as Wedzicha: cysteine + dicarbonyl capture is fast | not run |

Zhang 2025b's blank pot, in detail:
- TTCA plateau 17.7 / 18.2 mmol/L at 30 / 50 min; CAMA 1.61.
- Free Cys 0.80 and Xyl 1.22 at 50 min.
- ~19 mmol/L of each reactant sits in no measured pool. The deficit is equal for Cys and Xyl at every time and is
  12 - 14 mM already at 5 min (CF-063).
- Cys and Xyl print one identical ELSD curve.

### 2.5 Mass balance

| source | what it shows | for / against | engine |
|---|---|---|---|
| Zhang 2025b | 52 % of each reactant accounted for at 50 min | the gap | seen |
| Cao 2017; Zhao 2017b | over-close by 6 - 17 % and 9 - 16 %; cysteine-loss ceilings ~4e-4 and ~2e-4 /min at 90 C | against a large sink at 90 C | seen |
| Hou 2026 (Jiangnan; Xyl 0.2 M + Gly 0.1 M products with 5 mM Cys, 120 C step) | TTCA 0.193 mM (3.9 %) and cystine ~0.20 mM (~8 % as Cys equivalents) | cystine is a real but minor member (P3's band); most of the cysteine is in neither form | not run |
| Wei 2022 | 11 % of the consumed Cys is TTCA | as above | not run |
| Obata 1960; Mizutani 1960 (Hokkaido; 120 C bath, CO2 sweep, unbuffered) | free Cys gives H2S as 0.95 % of S in 6 h (k ~2.7e-5 /min); cystine lags, then exceeds Cys (x1.33 - 1.47 at 3 - 6 h); NH3 / H2S = 1.4 | the cysteine-only H2S sink is ~1e-4 of what Zhang's deficit needs | not run (Mizutani's F04#26 is E10's TEST contrast: item 7 only) |

### 2.6 One factor varied inside one study

| source | factor -> result | for / against | engine |
|---|---|---|---|
| Zhai 2021 Fig. 4b; Fig. 6 | label pH -> TTCA loss x10. Extra xylose speeds the loss modestly at pH 5.5 - 7 and slows it at pH 8 over the first 40 min | part B; release + recapture | seen |
| Kang 2026 SI Table S6 (storage, 60 days) | TTCA retention 95.51 / 96.70 / 88.81 % at pH 5.5 / 7 / 9 | non-monotone, <= 40 C: a storage cell never sets a processing law | not run; context |
| Zhang 2025b lysine arms | high-lysine free Cys is plotted below the ELSD curve's zero-signal point (0.223 mmol/L) | a low-end ELSD problem (R-cal) | not representable (no lysine in the lane) |
| Cao 2017 glycine arms | glycine speeds the TTCA drain | qualitative (3.2 - 9.6 mol/L glycine is not a solution) | not run |

### 2.7 Fitted time courses

| source | what it shows | engine |
|---|---|---|
| Schieberle 2000 Table IV (DFA; 100 C, 0.5 M phosphate pH 5) | FFT 2.0 / 3.1 / 110 / 132 ug at 30 / 60 / 360 / 720 min. R0 fits FFT almost exactly; C1 is ~2 dex low at 30 min. Signature (d) | **seen** (M1 D1) |
| B60's FFT rows (`thiol_submodel.yml`) | the FFT levels B60 was fitted on; C1 costs 36 nats on B60's 52 rows | seen |

### 2.8 Earlier engine records on the same parts

| record | part | result |
|---|---|---|
| supply screen (prereg v2 sec. 0B) | A, B, R-K, R-sink | S1, S2 and S3 fail (b) and (c). S3 + exit passes (a) and (b). No structure reaches (c) |
| M1, EX-061 | A | "not separated", artefact-driven. C1 "3 - 5x too slow" was read off Zhang's free-reactant decay, the data R-cal disputes, so it is **weak** motivation for freeing the release |
| B58, B59b | A | Cao's K by detailed balance; de Roos's route; on Zhai 2019 the CAMA formation must be 2 - 5x slower |
| CF-076 | A, R-K | no single pH-flat exchange holds dilute and concentrated pots. Kallen-shaped K(pH) admitted provisionally (P 0.25) |
| CF-063 | members, R-cal | the deficit is 1 : 1 and formed within 5 min; P(adduct rather than a destructive sink) ~0.65 |
| CF-107 | all | whole-model refusals against instant supply are not evidence against a reservoir |

### 2.9 The held-dossier screen (v2)

Every held dossier that mentions thiazolidine or TTCA (89 files; two reading agents, read-only, blinded) was
classified. Files not named elsewhere in this case:
- **Bearing, weakly or as context only:**
  - `du2024`: a dienal thiazolidine; its g/L values are flagged impossible.
  - `esterbauer1975`: an enal Michael adduct, 20 C.
  - `feng2022`: RG-ARP + Cys at 100 C; reacted Cys 9.15 of 20 mM, "data not shown".
  - `goldberg2002`: cysteine pKa and Delta H.
  - `haleva1999`: Cys + furanone, 70 C.
  - `hicks1975`: cites the TTCA pKa 6.24 only.
  - `k4c_wave_consolidated`, `k6a_sulfur_ladders_synthesis`: syntheses of the sources above.
  - `li2011`: Glc + Cys 0.2 M, free-Cys k 1.1e-3 - 2.0e-2 /min at 100 - 130 C; its assay may count the ring.
  - `li2016b`: an H2S method with no LOD.
  - `parker2014`: the reservoir proposed, not measured.
  - `wang2020`: GSH-sugar thiazolidines 0.8 %.
  - `yaghmur2005`: the furfural-Cys ring.
  - `yu2016`: ascorbate + Cys, 110 - 145 C.
  - `zhang2023`: Met / TTCA, volatiles only.
  - `zhang2026`: Cys + glucose -> Cys-Amadori, Ea 100.9 kJ/mol, already `r_cama_form`'s barrier; units unprinted.
  - `zhai2023jafc`: added Cys delays ARP and 3-DX, with no TTCA measured.
- **Not bearing:** `ames2001`, `bolton1993`, `chan1994`, `cui2022`, `davidek2006a`, `doorn2002`, `du2023`,
  `elmore1997`, `friedman1990`, `gao2020b`, `guan2012`, `hidalgo2004`, `hidalgo2016`, `hofmann1995b`, `kallen1973`
  (an erratum, no constant changes), `kerler2010`, `mandin1999`, `meng2017`, `mikami2017`, `nikolantonaki2012`,
  `qi2022`, `saunders1966`, `tai1997`, `tai1998`, `thornalley1999`, `vanseeventer2001`, `wang2026`,
  `weerawatanakorn2015`, `xu2008`, `xu2011`, `xu2013`, `yaghmur2002`, `zhang1989`, `zhou2024b`, `zyzak1995`. Each
  fails on one of these grounds: no cysteine or TTCA chemistry, end-point volatiles only, a matrix outside scope, a
  non-thiol sulfur source, a review, or an erratum.
- **Result:** no held dossier measures, at 90 - 140 C in water, any of the following:
  - a K for a cysteine-sugar thiazolidine other than Cao and Zhao 2017b;
  - a TTCA loss or release rate at held pH;
  - a Zhang-type cysteine and xylose closure by an independent assay.

  The cysteine-only sink measurements (Zhou 2022, Wei 2022, Li 2011, Liu 2022, Feng 2022, Yu 2016, Obata / Mizutani)
  bear on R-sink's cysteine arm, which cannot make a 1 : 1 deficit.

**Sources left out, with reasons.** These belong to other steps or regimes:
- the oxidant, O2 and disulfide parts (M1's C3, B11, B27);
- lipid, buffer and hexose cells (M4 / M5);
- the storage and low-water cells (context only).

Dossiers that mention neither thiazolidine nor TTCA were not screened.

## 3. The rivals

| id | explanation | (b) TTCA loss vs label pH | (c) Zhang's missing ~19 mM | (d) FFT shape |
|---|---|---|---|---|
| **R0** | instant supply, the incumbent (B60) | no TTCA in Cys + Xyl pots; fed TTCA's loss is pH-flat | no TTCA pool: reactants consumed in seconds | fits (the reference) |
| **C1** | M1's reservoir, release free, with K(T) banded as in every rival here (sec. 5.2); M1's T-flat C1 is reported beside it as the M1 reference | pH-flat (ratio <= 1.01 in the screen) | TTCA lost when free Cys is low; no pool | lost at 30 min unless a faster release restores it |
| **X** (the claim) | C1 + the exit (X-cama or X-deg) + reversible P1a, P1b, P2, P4a + irreversible P3, P4b | rises with pH; must also give the acceleration under the acid drift | what the members reach inside their bands: ~7 mM Cys and ~5 mM Xyl inside non-numerical bands (sec. 5.2) | open |
| **X0** | C1 + the exit, no pool | as X | not reached; Zhang's blank declared not representable | as X |
| **R-rel** | a base-favoured RELEASE (`k_ttca_cys` x 10^(n_r (pH_c - 5.5))) with K flat, so formation speeds up equally | TTCA re-forms unless free Cys is drained, so the loss rises only if the Cys sinks keep up | as C1 | open |
| **R-cal** | X0 + a **time- or handling-dependent analytical loss** of free Cys and Xyl after the quench, into an UNMEASURED product (a 1 : 1 loss into TTCA would raise the TTCA reading); mechanism unidentified, numerical band, which cannot be identified on Zhang because it prints no quench-to-injection times; variant R-cal+s adds a late sink | as X | no pool | as X |
| **R-sink** | X0 + material destroyed: S3's 1 : 1 drain on free Cys and PENT, plus a cysteine-only sink | as X | destroyed, never returned | open |
| **R-K** | Kallen-shaped K(pH) + the exit (CF-076's former lead) | K rising with pH stabilises TTCA, so the exit must outrun it | **fails by arithmetic**: the Kallen shape gives K(6.5) / K(4.4) = 3.4 - 18x for pKa_T 6.1 - 5.2, not the ~300x needed | open |

The columns need these notes:
- **R-sink** already failed (c) inside its band in the screen (S3: free Cys 15 / Xyl 18 mmol/L when TTCA and CAMA
  pass). It stays as the irreversible alternative for Arm 2, not as a candidate on disk.
- **R-K** is genuine only as CF-076's former lead.
- **R-rel** is against the analogues: Karpinski, Bustillo Trueba, Pesek and Nagasawa all show release not favoured,
  or stopped, in base.

**R-cal, the case for and against.**
- **For:**
  - At K = 0.060 L/mmol and TTCA 18.2 mM the free reactants should be sqrt(18.2 / 0.060) = 17.4 mM each. That closes
    the balance to (18.2 + 1.61 + 17.4) / 40 = 93 %.
  - Its strongest point: at full closure, Cao's K at Zhang's charge predicts a TTCA plateau of 18.3 - 20.1 mM for K
    0.045 - 0.060. That is the measured plateau (17.7 / 18.2). The K(T) band of sec. 5.2 widens this to 16.0 - 25.1
    mM (K(120) 0.032 - 0.144 L/mmol).
- **Against a static, level-dependent under-read:** it fails Zhang's own time course.
  - The 5-min deficit of 12 - 14 mM sits at free readings of 18 - 20 mM, where no low-end effect acts.
  - From 30 to 50 min, free Cys falls 4.0 -> 0.80 and Xyl 3.45 -> 1.21, while TTCA + CAMA rise only 0.5 mM.
  - None of the three named mechanisms gives an equal ~20x under-read of both reactants. A copied curve misreads one
    species. The 0.223 mmol/L zero-signal floor over-reads. A convex ELSD response under a linear fit misreads by tens
    of percent.

So R-cal is re-specified as a loss that grows with time or handling (between quench and injection, say), mechanism
unidentified, band numerical. The K-plateau agreement is evidence for it, and the time course is evidence against its
simple forms.

X, R-cal and R-sink can hold together in parts. The case asks which parts the decisive data need.

**The stop rule, restated (v3).** If no structure passes (b) and (d) within ~25 h, supply is "not identifiable on
disk". (a) is a consistency check (sec. 5.1) and does not enter the rule.

Feasibility and comparison are separate tests:
- **Feasibility (the stop rule's "pass") is a test on measurement noise and read error only, WITHOUT omega.** A
  structure passes a signature when both of these hold:
  - it meets the screen's declared tolerance (b: the 80-min ln-fraction ratio pH 8 / 5.5 >= half the measured 9.6;
    d: FFT 360 / 30 and 720 / 360 within the series' read-error band), or, equivalently, a posterior-predictive
    p > 0.05 under the measurement noise and read error alone;
  - it is not beaten on that signature by the pH-flat null (C1's exit, or R0's pH-flat `r_ttca_deg`) by more than the
    comparison threshold, 1.92 nats (the 90 % half-width on one degree of freedom).
- **Why not omega.** Under omega 0.419 a pH-flat rival scores p ~0.8 - 0.99 on (b), so C1 and even R0 would "pass"
  (b), and R0 passes (d). The rule could then never fire.
- **Which form governs.** The feasibility form governs the stop rule. Omega (sec. 5.1) is used only in the comparison
  between rivals.

On ranking versus deciding:
- **(b) can rank but not decide.** Against the pH-flat null it carries ~2.7 nats under omega, against the 1.92-nat
  threshold. That is a ranking, not a decision, consistent with the development-grade verdict.
- **(c) does not enter the rule.** It **ranks X, X0, R-cal and R-sink on development data only**; it cannot both
  decide and be unable to ship.

## 4. The decisive data

### 4.1 Why the screen's pots and the other seen data cannot decide

- **The screen's pots are seen.** The screen evaluated the engine on Cao 2017, Zhai 2021 and Zhang 2025b across every
  structure's band, and M1 scored Zhang's blank and Cao.
- **Cao also feeds constants.** It fixes K and anchors `r_cama_form`, and it set the screen's m and drain ceilings.
- **Other seen data:**
  - Zhao 2017b (M1's fit rows and side effects);
  - Zhai 2019 Fig. 3A / 3B;
  - Zhai 2023's Fig. 3a and Fig. 5;
  - Kang's Fig. S3 / S4 (and through them B60's `k_ttca_deg`);
  - Schieberle 2000 (D1) and B60's FFT rows;
  - de Roos 2005 Fig. 4 (spent as a constant).

A comparison on any of them can say only whether a structure is feasible, which is development grade.

### 4.2 Unseen data on disk: none is decisive

| source (not run) | what it could test | why it is not decisive |
|---|---|---|
| Hou 2017 (fed TTCA vs fed CAMA, 140 C, pH 5.5) | a lag in fed TTCA's thiols against fed CAMA's, same vial and detector (scale cancels) | one pH, so silent on part B and on the pool. Peak areas at 30 min only, and it needs keying. A side unit for X-cama vs X-deg at most |
| Wei 2022; Zhou 2022's Cys-alone arm; Li 2011; Liu 2022 | a cysteine-only sink | it cannot make a 1 : 1 deficit. Wei's sugar arrives through the GX-ARP's reversal, which the engine does not carry |
| Hou 2026 (TTCA and cystine after a 120 C step) | the cystine member's size | it sets P3's band |
| Kang 2026 end pH; Table S5; Table S6 | the pH trajectory; FFT / MFT vs pH; storage loss vs pH | the end pH is an input. S5 is compromised. S6 is storage |
| Zhai 2020; Zhai 2023c | TTCA -> CAMA vs pH (spray-drying); fed TTCA at pH 7.4 | outside the regime; no time course |
| Zhai 2019 Flag 1 ("TTCA disappeared completely" at 120 C with added Cys) | a reactive open-chain member | one qualitative sentence with no numbers, and it contradicts the same paper's Fig. 3A mass action. It names a member but does not size it |
| Wedzicha 1991; Nagasawa 1981; Radomski 1989; Kim 1986 | member and exchange classes | used as bands |

**No unseen decisive data exist on disk.** v1 put this at P ~0.85 before the full screen, and the review asked for
~0.75 until the screen was done. After the screen of sec. 2.9 it is **P ~0.85**. The residual is dossiers that never
name thiazolidine or TTCA but could hold a Zhang-type closure.

**The verdict of M1c is therefore development-grade.** A structure can be found feasible or infeasible on (b) and (d),
and the rivals ranked on (c), but nothing on disk can ship a supply structure as decided. The designed experiment in
4.4 would decide it.

### 4.3 The new papers (from the queue only, declared before reading)

| paper (queue description) | role, fixed now | the measurement that would make it decisive | which rivals it separates | P it holds that measurement |
|---|---|---|---|---|
| **belikov1967** (stability of L-cysteine + carbohydrate condensation products; possibly K for several sugar thiazolidines; title only) | **decisive only if at process temperature (>= 90 C)**; at 25 - 60 C a same-class band for K(pH) / K(T), by amendment | K (or ring : open ratio) for a cysteine-PENTOSE ring at >= 2 pH or >= 2 temperatures, >= 90 C | flat K (C1, X, R-cal) vs a pH-dependent K (R-K; Zhao 2017b's drift); P1a's T-law | ~0.08 |
| **wlodek1993** (non-enzymatic hydrolysis of sugar TCAs; a bis-Cys thioaminal pool) | **decisive only if at >= 90 C**; otherwise (the queue implies 37 C) a band source for R-rel's release pH law | a pentose TCA's release rate at >= 2 pH between 6 and 8, cysteine measured | R-rel vs X (is the base-favoured loss a release?). A bis-Cys pool (2 Cys : 1 aldehyde) cannot be Zhang's 1 : 1 deficit | ~0.15 |
| **davidek2003** (HPAEC-ECD / DAD of sugar + amino acid + Amadori together; detection floors) | **decisive for R-cal only on a Zhang-like pot** (cysteine + pentose, ~40 mM, pH ~6.5, ~120 C); otherwise floors for its own assay | free Cys, free pentose and TTCA / CAMA by an assay independent of ELSD, >= 2 times | R-cal vs X / R-sink | ~0.05 |
| gajda1990 (protonation constants of 2-(polyhydroxyalkyl)TCAs, pH 1.5 - 8) | band source: TTCA's ring-N pKa (R-K) | -- | -- | -- |
| canle1996 (N-aryl thiazolidines, ring-opening k and ring / open K vs pH, 25 C) | band source, distant class (N-aryl): the shape of P1a's K and of ring opening | -- | -- | -- |
| szilagyi1979 (C-2 inversion and C-4 epimerisation in TCAs) | band source: the P1 / P2 exchange speed; the epimer share | -- | -- | -- |
| davis1991, davis1991b (penicillin thiazolidine ring opening, pH-rate profiles) | band source: the pH shape of ring opening (R-rel) | -- | -- | -- |
| carta1996 (L-cystine solubility vs pH and NaCl) | band source: P3's ceiling; whether cystine precipitates on cooling (an operator question) | -- | -- | -- |
| sohn1995 (ammonia from amino acids heated in water) | band source: R-sink's cysteine-only arm; CF-084 (M1b) | -- | -- | -- |
| davidek2004 (pentose-Amadori assay; xylose + glycine) | evidence only (glycine, not cysteine) | -- | -- | -- |

The combined chance that at least one of the three candidates supplies a decisive measurement is **~0.25**
(1 - 0.92 x 0.85 x 0.95 = 0.26). None of the eleven sets a band in this case. A band source enters only through an
amendment registered before any engine result of M1c exists.

### 4.4 The designed experiment (W6)

The design is one laboratory, three arms, with these conditions throughout:
- n = 3 independent heatings per cell;
- sealed vessels with stated fill and headspace (N2 and air arms in Arm 2);
- pH measured on the cooled scale for every sample, and in situ where an electrode allows;
- every assay's LOD / LOQ and spike recovery reported;
- the time from quench to injection logged for every sample.

| arm | conditions | species | time points | decides |
|---|---|---|---|---|
| **1. the exit's pH law at held pH** | purified TTCA (equilibrium epimer mixture, NMR-characterised) 10 mmol/L. pH 4.5, 5.5, 6.5, 7.5, 8.0, held by pH-stat or in phosphate at 0.05 and 0.2 M (extrapolate the buffer term to zero). 100 and 120 C (90 C at pH 4.5 bridges to Cao). +- 10 mmol/L [13C5]xylose | TTCA (each epimer); CAMA (each isomer); free Cys derivatised at the quench (N-ethylmaleimide, labelled internal standard); free xylose (HPAEC-PAD); [13C5]-TTCA; cystine; H2S (trap); 1-DX / 3-DX; A420; pH. **Quench control:** TTCA standards through the same NEM quench at 0 / 10 / 30 / 60 min before injection, to show that the trap does not drain TTCA through its open chain (Schubert's iodoacetate trap freed cysteine from a thiazolidine at room temperature within ~0.5 h) | 0, 2, 5, 10, 20, 40, 60, 90, 120 min | n at fixed pH, with no drift. Which step carries it: CAMA appears (X-cama); deoxyosones or browning without CAMA (X-deg); free Cys and xylose rise transiently, with the xylose recovered (R-rel; the cysteine sinks drain the Cys, so it does not stay). An induction at constant pH means a P1-type pool |
| **2. Zhang's pool in a second laboratory** | Cys + Xyl 40 + 40 mmol/L, 0.2 M phosphate pH 6.5, 120 C, 5 mL in 15 mL (Zhang's design). The 0.4 M / 90 C recipe (Cao's) as a bridge. A [13C5]xylose run | free Cys (as in Arm 1); free xylose (HPAEC-PAD); **split samples run by Zhang's own HPLC-ELSD method, with separate matrix-matched Cys and Xyl calibrations in 0.2 M phosphate over 0.2 - 40 mM**; TTCA (**epimers resolved**) and CAMA by LC-MS/MS with labelled standards; total releasable cysteine (a validated exchange checked on TTCA standards); cystine in solution and in any pellet; 13C NMR of [13C5]xylose samples quenched at 10 and 50 min; the 50-min sample diluted 20x at 120 C and followed for 30 min; cold dilution vs 4 C storage at 0 / 4 / 24 h | 0, 1, 2, 5, 10, 20, 30, 50 min | (c), as in the decision table below |
| **3. K(pH) and K(T)** | Cys + Xyl equilibrated from both sides (from TTCA, and from Cys + Xyl) at pH 4.5, 5.5, 6.5, 7.5 and 25, 60, 90 C. At pH >= 6.5 above 60 C the exit drains TTCA during equilibration: there K is fitted with the exit in the model, or high-pH K is limited to <= 60 C. TTCA ring-N pKa titrated at 25 and 90 C | all four species, and pH at temperature | to a plateau (two successive samples agree) | flat K vs R-K and Zhao's drift; P1a's T-law; K(T) |

**The Arm 2 decision table** (directions from mechanism):

| measurement | X, reversible members | X, irreversible members | R-sink | R-cal |
|---|---|---|---|---|
| closure by the independent assays (HPAEC, labelled LC-MS/MS) | < 100 % (the pool is unmeasured), closed by NMR or by total releasable Cys | < 100 %, not closed by releasable Cys | < 100 % | **~100 %** |
| Zhang's own ELSD method on split samples | reproduces the 1 : 1 deficit | reproduces it | reproduces it | reproduces it **only** where the independent assays close: the loss is in that method |
| 20x dilution of the 50-min sample | **Cys + Xyl returned beyond TTCA's own dissociation** at K | only TTCA's own dissociation | only TTCA's own dissociation | only TTCA's own dissociation |
| 13C NMR at 10 / 50 min | open-chain, xylosylamine or epimer signals of the size of the gap | adduct or cystine signals | no pool signals | no pool signals |
| quench-to-injection time | no trend | no trend | no trend | **a deficit that grows with the delay** (the time-dependent form) |

The arm can refute R-cal (independent closure fails) but confirm it only jointly with a delay trend or a
method-specific deficit.

Request the experiment through W6, with CF-076's designed experiment (TTCA 10 vs 100 mM, +-phosphate, +-morpholine)
folded into Arm 1.

## 5. Tolerances from the data

### 5.1 The measurement model per unit

**Rules (M1.1).**
- tau (H0's tau_ABS 0.964 dex) enters only on absolute levels compared across laboratories or against an engine
  level.
- omega = **0.419 dex = 0.965 ln** in every unit's primary covariance, per printed value. This is H0's ship value: of
  the two sets (reset 0.390, warm 0.419), the worse one governs. 0.390 is a sensitivity.
- Floors come from each assay's LOD / LOQ or calibration range, never from a data point.
- Densities are exact or Monte-Carlo (PSIS k-hat < 0.7, ESS >= 100) on every form that can flip a sign.

| unit | role | value basis | reported noise | error terms in the primary | floor | note |
|---|---|---|---|---|---|---|
| (a) Cao 2017 and Zhao 2017b, log K_app vs measured pH | **consistency check only; not a likelihood unit** | a within-sample quotient of three absolute levels | n = 3, RSD < 10 % per species; read errors +-0.004 / 0.005 mol/L; +-0.03 pH | omega once per quotient: slope sd **0.95 dex / pH** (0.89 at 0.390) over Cao's pH spread | none needed | it separates no rival. C1, X, X0, R-cal and R-sink are flat by construction at pH ~4.4, and R-K moves <= 0.075 dex there. Reported: whether each structure's K_app stays inside Cao's and Zhao's measured ranges (58.5 - 64.6; 59 - 78 L/mol). Zhao's -0.29 dex / pH sits on a 6-fold-dilution pH scale |
| (b) Zhai 2021 Fig. 4b | stop-rule unit | **primary: matched-time difference ln(C_pH(t) / C_5.5(t))** for pH 6, 7, 8 and the +Xyl series against their own pH 5.5. Any scale within one assay cancels; C0 is likely the nominal 10.00 mM | SD (n = 3) for TTCA alone at pH 5.5 / 6 / 7 only. For pH 8 and the +Xyl series, a transfer: the median printed relative SD of the same figure (~3 - 4 %), sensitivity x0.5 / x2. Read error +-0.005 mmol/L | omega per printed value (2 omega^2 per difference). **Sensitivity: a shared within-series error term** (omega split into a series-level part and a point part, half each), since independent omega per point overstates power | TTCA LOD 2.21 mg/L = 0.0087 mmol/L (Zhai's printed); "0.00" points censored | arithmetic on the measured series: against a pH-flat rival the pH 7 and 8 differences carry chi^2 5.5, **~2.7 nats** under omega 0.419. Modest, and less with the shared term: against the 1.92-nat threshold (b) can rank but cannot decide. Its stop-rule feasibility form is on measurement noise only (sec. 3). The pH 7 points at 100 / 120 min print SD 0 (not measured as SD): the transfer applies |
| (b) product routing, Zhai 2021 Fig. 5 | **declared feasibility check** (not a likelihood unit) | standing 3-DX and 1-DX at pH 8, 100 C, 120 min, semi-quantitative (diacetyl surrogate) | as printed | -- | -- | a variant passes if its predicted 3-DX is <= 3x of 0.8 % of the TTCA lost and its 1-DX <= 3x of 106 umol/L. The factor 3 is the surrogate's declared response-factor uncertainty, left as declared. A failure is reported; the check vetoes nothing |
| (c) Zhang 2025b | ranking unit (development) | **primary: the closure fraction** (Cys + TTCA + CAMA, and Xyl + TTCA + CAMA, over the 40 mM charge) at each time | n = 3, error bars small and not transcribed; read errors +-0.3 (TTCA), +-0.05 (CAMA) mmol/L | omega per printed value, plus a **per-assay calibration term**: TTCA 8 %, CAMA 20 %, Cys / Xyl ELSD 15 %. These are the Cao-vs-Zhao 2017b same-lab, between-run and between-assay spreads (zhao2017b dossier): noise terms, not model bands | **Zhang's own: 0.223 mmol/L** (the printed ELSD curve's zero-signal point) for Cys / Xyl; **0.1 mmol/L** for TTCA and CAMA (the lowest standard). Zhai 2021's 0.025 / 0.020 (UNVERIFIED for Cys / Xyl) is a sensitivity only | the K gap: Zhang's K_app 18.6 L/mmol against Cao's 0.060 is 2.49 dex. Data against data, on measurement noise only (~0.1 dex), that is ~25 sd. A flat-K model against Zhang's K_app needs one quotient's model error: **3.4 sd** with omega per printed value, **5.9 sd** with omega once (3.7 / 6.4 at 0.390). A large, one-pot discrepancy |
| (d) Schieberle 2000 FFT and MFT ratios 360 / 30, 720 / 360 | stop-rule unit | ratios within one series (no tau); bivariate, sharing the 360-min point | no SD or n printed; M1's 0.094 ln per point | omega per point | 1 ug/L | as M1 and prereg v2 |
| (d) B60's FFT rows | in-sample guard | absolute levels, several labs | each row's printed SD or transfer | omega + tau | each row's floor | -- |
| decisive, Arm 1 | decides | ln(C(t) / C0) per species and pH; n from the within-study pH contrast | the lab's replicate SD | omega per point; a shared within-series term as sensitivity | the arm's LOQs, censored likelihood | exact densities |
| decisive, Arm 2 | decides | closure fraction per time; the dilution return ln(C after / C before); the delay trend | replicate SD; spike recoveries | omega per printed value; a per-assay calibration term from the arm's own recoveries; tau only against Zhang's own levels | the arm's LOQs | isotope dilution keeps the calibration term small |
| decisive, belikov1967 (at >= 90 C) | decides | ln K ratio between two pH (or T) values in one study | as printed; else a declared transfer (RSD 10 % per species), sensitivity x0.5 / x2 | 2 omega^2 on the ratio (omega once per quotient) | the paper's floors | -- |
| decisive, wlodek1993 (at >= 90 C) | decides | ln ratio of release rates between pH values in one study | as printed; else the transfer above | omega on the underlying levels | the paper's floors | -- |
| decisive, davidek2003 (Zhang-like pot) | decides | closure fraction in one study | its printed RSDs and recoveries | omega per value; a per-assay term from its recoveries | its printed LODs | -- |

### 5.2 The bands and their sources

No band below comes from Cao 2017, Zhao 2017b, Zhai 2021, Zhang 2025b or Kang 2026 (all screened or seen), or from any
candidate named in 4.3 (decisive). Where the only source is a seen pot, the value enters as a fixed measured constant
or as a required-agreement sensitivity, never as a band, and the case says which unit it then cannot test.

| coordinate | band | source | kind |
|---|---|---|---|
| n, the exit's pH slope (dex per pH unit, instantaneous cooled pH) | **0.30 - 0.74** | re-derived from Martins & van Boekel 2003 Part II Table 3 (see the note below this table) | same-class transfer; the log-linear form is an assumption (sensitivity 0 - 1: one hydroxide) |
| X-cama level | `r_cama_form` 3.29e-3 /min at 90 C, Ea 100.9 kJ/mol, referenced at pH_ref 4.4 | B59b anchor on Cao's CAMA rise (seen); Zhang 2026's Cys + glucose barrier (hexose transfer). Zhao 2017b's floor ~7e-4 /min is a required-agreement sensitivity (the anchor x0.2) | measured constant from a seen pot: Cao then tests only K's consistency, not the CAMA rise |
| X-deg level | `k_ttca_deg` at its B60 value +- 1 dex, pH_ref 5.5 | the centre is B60's fitted value, which rests on Kang's 16.3 % release ceiling (= Zhai 2021 Fig. 6), a seen pot; the width has no source | numerical; the centre is declared as from a seen pot |
| `k_ttca_cys`, the release | log10 k(145 C) = B60 +- 3 dex; Cao's floor (>= 7e-3 /min at 90 C) as a required-agreement sensitivity only | stage-1 convention; Cao arithmetic 4 | numerical (primary) |
| R-rel: n_r, the release's pH slope | 0 - 1 | mechanism (one hydroxide); the analogues (Karpinski, Bustillo Trueba, Nagasawa) all say <= 0 above pH 6 | numerical; the analogues argue against any positive value |
| K at 90 C | 0.060 L/mmol fixed; sensitivity **0.045 - 0.080** | Cao 2017 and Zhao 2017b: the only two cysteine-pentose K records on disk, from one laboratory and one recipe (Zhao's 59 - 78 widens the upper end) | measured constants from seen pots |
| K(T): Delta H of ring formation | **-25 to +35 kJ/mol**, giving K(120) / K(90) 0.53 - 2.4 | Kim 1986 (aliphatic thiazolidine, +4.9 kJ/mol, 13 - 38 C) +- 30 | same-class anchor, numerical width. Kim's Delta H is measured against a FREE aliphatic aldehyde. A sugar's apparent K also carries the sugar's (positive) ring-opening enthalpy, so the band's upper end may be low. The R-cal arithmetic and the members' sizing assume K(120) = K(90) and are re-stated over this band |
| R-K: pKa_T | 5.2 - 6.1 | Radomski 1989 (inference, as CF-076); gajda1990 replaces it by amendment if read | inferred class |
| P1a: f_oc = [open] / [TTCA] | 0 - 0.05 at <= 37 C; Delta H_open 0 - 60 kJ/mol | Radomski 1989 (not seen by NMR, 37 C); Karpinski 1990 (no open thiol, LOD ~1 uM) | same-class at <= 37 C; **numerical** in T |
| P1a / P2 exchange speed | >= 2e-3 /min (neutral, room temperature), acid-catalysed | Nagasawa 1981 (penicillamine analogue); Radomski 1989 (~0.2 /min, 37 C) | class bound: fast on the pots' timescale |
| P1b: N-xylosylamine share | 0 - 0.05 of TTCA | none measured; set equal to P1a's (both are open or non-thiazolidine 1 : 1 forms) | **numerical** |
| P2: epimer share not integrated | 0 - 0.20 of TTCA | Radomski 1989 (2R : 2S = 4.06, 37 C) | measured (TTCA itself); operator switch |
| P3: Cys -> 1/2 cystine | first-order ceiling 7e-4 /min at 120 C (<= 8 % in 120 min); barrier 40 - 135 kJ/mol | Hou 2026 (Jiangnan, 120 C step); `sugar_trunk` class | same-network ceiling; numerical barrier |
| P4a: second CAMA isomer outside the MRM | 0 - 1 of CAMA | Zhai 2021 SI and Zhai 2019 NMR show two CAMA forms; no share measured | numerical; a verdict that needs it does not ship |
| P4b: 3-DX + Cys (`r_tdp_cys_trap`) | binding fast, 1 : 1 (not rate-limiting): ~0.016 /min at 10 C if first order in the osone (untested; zero order in Cys over 10 - 50 mM). Size capped by the 3-deoxypentosone supply | Wedzicha 1991 (3-deoxyhexosulose, pH 5.5); Edwards 1992 | same-class (hexosulose for pentosulose) |
| R-cal: the time- or handling-dependent loss | 0 - 0.1 /min of the measured free Cys and Xyl after the quench, into an unmeasured product | none on disk (hu2012, queued and unread, is a candidate class source for ELSD non-linearity). Not identifiable on Zhang, which prints no quench-to-injection times | **numerical** |
| R-sink: the 1 : 1 drain | 1e-6 - 3e-3 /min at 90 C; barrier 40 - 135 kJ/mol | numerical; Zhao 2017b's ceiling 2e-4 /min at 90 C as a required-agreement check | numerical |
| R-sink: the cysteine-only sink | 0 - 2e-4 /min at 90 C; at 120 C, Zhou 2022's ~2.5e-3 /min (Cys alone, water) as the upper corner | Zhao 2017b's ceiling (seen): sensitivity; Zhou 2022's Cys-alone arm (not run); sohn1995 by amendment | numerical primary |
| pH trajectory, Zhai-design pots (cooled scale) | **three declared shapes, each run:** (i) no drift (label pH held: the bound that favours a base-favoured exit most); (ii) linear in time from the label to Kang's end pH at 120 min (5.1 / 4.9 / 4.7; 5.0 interpolated for 6); (iii) front-loaded, an exponential approach to Kang's end pH with tau 20 min | Kang 2026 p. 3242 for the end points (its temperature unresolved). The front-loaded time scale is from Hou 2026's measured unbuffered trajectory (Jiangnan, 2 pH units in 20 min, a different pot). The linear shape has no source | measured end points; declared shapes. A verdict that holds under only one shape is reported as shape-conditional |
| pH trajectory, Cao 2017 | the keyed cooled readings (Fig. 4c) | Cao 2017 | measured |
| pH, Zhang 2025b | label 6.5 (0.2 M phosphate; no reading printed); sensitivity: the engine's charge-balance pH | -- | declared |

**Note on n's band (re-derived from Martins & van Boekel 2003 Part II Table 3).**
- **The steps:** DFG (glucose-glycine Amadori) decay steps k1 - k3 (`trunk_conditions.PH_SOURCE`).
- **The re-derivation:** the ratios k(6.8) / k(5.5) are 3.0 / 15.6 / 8.6 at 100 C and 8.0 / 7.3 / 9.8 at 120 C, on
  the INITIAL label pH of drifting pots.
  - Measured drops: 0.8 at 100 C and 0.9 at 120 C from pH 5.5; 0.19 and 0.21 from pH 6.8 (martins2003 dossier sec. 2
    and flag 7).
  - Under linear drift, the time-averaged gap is 1.60 (100 C) and 1.65 (120 C) units, not 1.3.
  - That gives 0.30 / 0.74 / 0.58 at 100 C and 0.55 / 0.52 / 0.60 at 120 C.
- **The window:** pH 5.5 - 6.8, extrapolated here to 4.15 - 8. TTCA's ring-N pKa (5.2 - 6.1) lies inside that range,
  so a log-linear pH law is an assumption, not a transfer.
- **Correction:** v1 used the 0.37 - 0.92 of `trunk_conditions.py` (1.3-unit gap).

**What the bands allow (arithmetic, no engine), at K(120) = K(90).** Inside the non-numerical bands, at Zhang's
plateau:
- P1a holds at most ~0.9 mM (5 % of 18 mM).
- P2 holds at most ~4.5 mM.
- P3 holds at most ~1.4 mM of cysteine by 50 min (7e-4 /min).
- P4b holds at most ~1 mM, capped by 3-deoxypentosone supply: standing 3-DX in comparable pots is <= ~1 mM (Zhai
  2021 Fig. 5; Hou 2026), and 12 - 14 mM in 5 min is out of reach.
- The total, ~8 mM of cysteine and ~6 mM of xylose, falls short of the 19 mM of each that (c) needs.

Over the K(T) band the shortfall stays. So X reaches (c) only through a numerical band (P1a at high T, P1b, P4a), or
not at all: P(X passes (c) inside non-numerical bands) ~0.1. R-cal's form needs a numerical band too. (c) ranks; Arm 2
decides.

## Open questions for the build

1. **The cooled scale.** n is banded on cooled-scale pH (Martins, Cao and Kang report that scale), but the chemistry
   acts at the in-situ pH. The primary uses the engine's cooled-scale pH of the current composition
   (`observations.ph_cooled`); the in-situ argument is a required-agreement sensitivity.
2. **X-cama against de Roos and Zhai 2020.** Both variants run; neither source chooses.
3. **The acceleration.** If no variant produces Zhai 2021's acceleration on any of the three pH shapes, (b) is still
   judged on the matched-time differences; the shape is reported.
4. **Zhao 2017b's rising K.** It is on a 6-fold-dilution pH scale and a different assay. If both records are read as
   replicates of one pot, K_app is flat within ~0.1 dex to 2 h and diverges after (TTCA by ELSD vs LC-MS). Is the
   divergence assay or chemistry? Only Arm 3 decides.
5. **Keying Hou 2017** (2 - 3 h, two readers) as a side unit for X-cama vs X-deg before the screen, or leave it out
   with the reason above.
6. **The shared within-series term** in (b): its split (half series, half point) is declared, not measured.

## Revision log

- **v1 (88457d7e), 2026-11-05.**
- **v2, 2026-11-05**, after an independent blind review (FREEZE AFTER FIXES). Still before any result. No new paper
  was opened, and no engine was run.
  - Zhao 2017b added (K drifts 59 -> 78 L/mol; recomputed slope -0.290 dex / pH, 90 % CI -0.331 to -0.248); the K
    sensitivity widened to 0.045 - 0.080; Zhao's ceilings used for R-sink and X-cama.
  - Added Wedzicha 1991 (P4b's band), Nagasawa 1981, Lienhard 1966 (excluded by arithmetic), Schubert 1936 / 1939,
    Roberts 1987 and Kim 1986 (no usable constants, except Kim's Delta H for K(T)), Obata / Mizutani 1960, Zhang
    2023b, Edwards 1992, Wondrak 2002, Pripis-Nicolau 2000, Zhou 2022 and Liu 2022.
  - The 89-dossier screen (sec. 2.9).
  - R-cal re-specified as a time- or handling-dependent loss, with its K-plateau agreement and its time-course
    failures stated.
  - R-rel added. The N-xylosylamine (P1b) added. A K(T) band added.
  - Members split into reversible and irreversible.
  - R-K's arithmetic failure (3.4 - 18x) and S3's screen failure stated.
  - n re-derived on the time-averaged pH gap (0.30 - 0.74).
  - X-cama referenced at pH 4.4. X-deg's seen centre declared.
  - (b)'s primary is now the matched-time difference, with a shared-error sensitivity and the routing check. (c)'s
    primary is now the closure fraction with per-assay terms, and Zhang's own floors.
  - (a) demoted to a consistency check. The stop rule restated on (b) and (d), with "pass" defined.
  - The K-gap wording corrected.
  - The designed experiment gains split ELSD samples, matrix-matched calibrations, quench-to-injection times,
    resolved epimers, a NEM quench control, the high-pH K restriction, the dilution prediction and the decision
    table.
  - omega 0.419.
  - New-paper roles conditioned on temperature and pot; combined P 0.3.
  - P(no unseen decisive data) re-stated after the screen.
- **v3 (FROZEN), 2026-11-05**, after the confirmation pass. Still before any result. No new paper was opened, and no
  engine was run.
  - N1: feasibility (the stop rule's "pass") is now on measurement noise and read error only, without omega, and a
    passing structure must not be beaten by the pH-flat null by more than 1.92 nats. Omega is kept for the comparison
    between rivals. (b) ranks but does not decide.
  - N2: Kim's Delta H caveat added; C1 carries K(T) banded; R-cal's lower plateau corrected to 16.0 mM.
  - N3: R-rel's Arm-1 signature corrected.
  - N6: Zhou 2022's Cys-alone arm marked not run.
  - N7: R-cal's loss goes to an unmeasured product, and its band is unidentifiable on Zhang; the routing check vetoes
    nothing.
  - m8: combined P ~0.25.
