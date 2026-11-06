# M0 integrity batch 1 (Strecker, sulfur): re-reads and the CK-01 / 06 / 07 / 02 checks

> **Misfit protocol step 1 (integrity)** ([OPERATING_PROCEDURE.md](../../docs/guides/OPERATING_PROCEDURE.md) sec. 4)
> for the Strecker and sulfur misfits of the misfit inventory (PR #82, branch `misfit-inventory`,
> `docs/reference/misfit_inventory.md`, batch 0 and sec. 6). 2026-11-06, from `origin/main` at `4174ed8d`.
> One object per check in [m0_integrity_batch1_strecker_sulfur.json](m0_integrity_batch1_strecker_sulfur.json).
>
> **Rules kept.** Each check's pass rule was written before the check ran, and every outcome is recorded, including
> "integrity holds". Sources were re-read from the PDFs in the main checkout (laptop only; cited here by dossier and
> page, never by path). No fit; nothing written under `data/`, `src/`, the registers, cases, preregs or the inventory.
> Engine use was read-only: `predict()` on the shipped engine, 1 core, niced.
>
> **Blinding.** Not read: hofmann2000b Figs. 2-3 and the Results text that discusses them (the Methods paragraph of
> the 85 C trapped experiment was read as text only, for its atmosphere wording); the B59b hold-outs; the E7 / E13 /
> E14 TEST sets; the WUR samples.

## Outcomes

| check | misfit | outcome |
|---|---|---|
| CK-07 | MF-17 | **No governance bug.** The Zhang 2025b ACTZ row is a declared blind hold-out, scored once as declared; it is now spent |
| CK-01 | MF-29, MF-33 | **Integrity holds.** All 12 Table 2 cells, the basis and every condition statement agree with the page |
| CK-50 | MF-33, MF-29 (d) | **Inconclusive** (headspace unprinted). R2 needs headspace / liquid r >= 2.9 (air arm), >= 5.8 (air + Cu arm, which is air) |
| CK-06 | MF-28 | **pH is passed; the misfit is structural**, and located: alanine routes the pot to the acrylamide lane, which has no pH law on the trunk steps |
| CK-02 | MF-21 | **Integrity holds; not a zero prediction.** A within-paper MGO ratio error (~1.27 dex) seen through a shared free scale |
| re-read | MF-10 | **Integrity holds.** Givaudan ribose / NF MFT = 43.5 / 0.9 = 48.3 per mol Cys |
| re-read (CK-10) | MF-15 | **Integrity holds; not explained.** The charges are right and reach the engine; the thiamine route has no cysteine term |
| re-read | MF-16 | **Integrity holds.** The share is molar on both sides; one 0.004 dex target slip; a calibration-intercept lead (step 3) |
| re-read | MF-30 | **Integrity holds.** A heat-up lag would widen, not explain, the early deficit |
| re-read | MF-33 | **Integrity holds** on the source; B29's engine pots are declared glycine stand-ins |

No misfit is explained by integrity in this batch. No keying error was found in a curated file under `data/`.

---

## CK-07: the Zhang 2025b ACTZ row (MF-17)

**Pass rule (before the check).** A governance bug if the row (a) enters any fit objective or calibration (fit rows,
the fit-target ledger, the B60 or R3En objective, a shipped constant), (b) was scored before Amendment 45 declared it
blind, or (c) sits in a contrast-split TEST partition while being used for development. Scoring it on the panel's
hold-out section after Amendment 45 is the hold-out's declared test, not a breach.

**Outcome: no governance bug.**

- The bundle `mp_holdout_xylose_cysteine_phosphate_120C_50min_pH65_Zhang2025b` was added in `0e6e6097`
  (2026-10-28 01:21, Amendment 45). It was first scored in `ab1fcc55` (2026-10-28 12:31, "B60 completed ... scorecards
  on B60"). So the bundle came before the score, as the amendment requires ("the first scorecard that reads these
  rows is their test").
- `core_panel_scores.md` scores it in panel `maillard_path_holdout`, tier REFERENCE, class `external_holdout`
  (ACTZ 114.13 vs 0.065 ppb, 3.245 dex). That is the hold-out section, not an objective.
- No fit reads it. `k_cys_actz`'s offset (R3En -2.64) is fitted on Zhou 2023's pH 7 row only
  (`fit_joint_deoxyosone_thiol.py`: `ACTZ_TABLE = thiol_zhou2023_actz`, `ACTZ_FIT_PH = 7.0`). No zhang2025b ACTZ row
  is in `data/lit/fit_rows`, `fit_targets.py` or `kinetic_core_b60_fit_targets.json`.
- The split files (`split.yml`, `split_v2_per_change.yml`) carry zhang2025b only as F04#19 (TTCA paired arms with
  lysine, DEVELOPMENT). They have no ACTZ cell, and zhang2025b is not a v2 calibration input.
- CF-009 (Table 2 footnote: "sums detected from both the DB-WAX and HP-5 MS capillary columns") can make each level
  at most 2x (0.30 dex) high. It cannot carry a -3.25 dex miss.

**Governance notes (no fix made here).**
1. The bundle's `fit_holdout_status` ("INDEPENDENT ... checked by grep on 2026-10-28") is stale. Since then the
   same pot's intermediates (Fig. 2 TTCA, CAMA, Cys, Xyl) have been read by:
   - `route_submodel.yml` (an unshipped sub-model);
   - B58 T1, the B59 hold-out diagnosis, and B60 verification V2 / V3;
   - the supply signature screen, and the M1 and M1 + M2 preregs.

   None of them ships a constant, so the ACTZ row is still independent of shipped constants. The coordinator should
   reword the status to "row-independent; same-pot intermediates seen at design" (a `data/` edit).
2. The inventory's MF-17 line "hold-out status is inconsistent ... yet the ACTZ row is on the panel" should read:
   "a declared blind hold-out, scored once as declared (2026-10-28); now SPENT. No ACTZ route candidate (R37-R39,
   the NH3 thiazole route) may claim it as out-of-sample confirmation."

## CK-01: Hofmann & Schieberle 2000b Table 2 and Methods (MF-29, MF-33)

**Pass rule (inventory).** Integrity fails if any keyed cell, basis or condition differs from the page. If the argon
cell's basis differs, MF-29 is explained if the corrected datum lies within the rivals' q10-q90.

**Read:** the text layer of pp. 4301-4302 and a render of p. 4302 (Methods and Tables 1-3 only; text and render
agree). Compared with `strecker_hofmann2000b.yml` and the dossier `hofmann2000b_extraction.md` (secs. 0-3, Table 2
and 3 parts).

| item | page | keyed file / dossier | agrees |
|---|---|---|---|
| Table 2, ARP-Phe PA / PAA (umol/mmol) | argon 0.6 / 0.3; air 5.5 / 3.0; air + Cu 13.8 / 7.6 | rows = printed / 10 mmol/L | yes |
| Table 2, Glc/Phe PA / PAA | 0.4 / 0.2; 1.4 / 1.8; 2.6 / 3.9 | same | yes |
| basis | title "Concentrations (umol/mmol)"; footnote a "precursors (1 mmol each) ... (10 mL ...)" | 1 umol/mmol = 0.1 mmol/L; argon PA + PAA 0.09 mmol/L | yes |
| scale | Experimental 0.1 mmol in 1 mL; footnotes 1 mmol in 10 mL | both declared; 100 mmol/L either way | yes |
| phosphate | Experimental "0.5 mmol/L"; footnotes (Tables 1-3) "0.5 mol/L" | 0.5 mol/L taken, CF-020 | yes |
| vial | "closed vial"; size, fill and headspace not printed | "size and headspace not stated" | yes |
| deoxygenation | not printed; "under an atmosphere of argon" (fn a), "pure argon" (text) | "not stated" | yes |
| air arm | fn b "Argon was replaced by air oxygen"; text "Flushing of the models with oxygen (expt II)"; fn c "air atmosphere" + 0.05 mmol CuSO4 | air taken, pure O2 the alternative | yes |
| ARP-Phe purity | "purity of about 90%" (24 mmol, 8 % yield) | declared (alternative 90 mmol/L) | yes |
| quantification | SIDA, [13C2]PA and [13C2]PAA, 5 mL water, pH 3.0, ether, GC-MS/CI | SIDA | yes |

**Outcome: integrity holds.** MF-29 is not explained by integrity.

**One sharpening of the air reading.** The Methods use the same phrase for the 85 C experiment: "under an atmosphere
of argon or air oxygen in the presence of copper(II) ions". The analogous Table 2 arm (footnote c) says "air
atmosphere". So "air oxygen" is the authors' word for air (German "Luftsauerstoff"), and "flushing ... with oxygen"
reads as flushing with air. P(the air arm is air) ~0.85 (mine).

**Proposed CF (MF-29 has none).** `kind: open_condition`, `sources: [hofmann2000b]`:

> Hofmann & Schieberle 2000b Table 2, expt I (fed ARP-Phe, argon): PA + PAA 0.9 umol/mmol (0.09 mmol/L) as printed.
> Re-read 2026-11-06 (`m0_integrity_batch1_strecker_sulfur`): the cells, the basis and the conditions agree with the
> page. Every M3 rival predicts 0.39-0.48 mmol/L (+0.63 to +0.73 dex), and integrity is ruled out.
>
> Three conditions that bear on the argon arm are not printed:
> - the deoxygenation method. Residual O2 would RAISE the argon datum, so it cannot explain an over-prediction;
> - the vial headspace;
> - the reaction scale (1 mL in the Experimental, 10 mL in the footnotes).
>
> Open explanations: the anoxic exits feeding route A (e1 / e3 / eF) are too fast; k_d (CF-121); s_DG1; ARP purity
> (<= 0.05 dex). Decider: DM:strecker_routes (a), with the headspace and the deoxygenation stated.

## CK-50: the O2 arithmetic (MF-33 (a), MF-29 (d))

**Pass rule (inventory).** The unit is umol O2; the error is the rounding of the printed vial volume. Under air, R2 is
refuted if the available O2 is < 250 umol. If the headspace is unprinted, the result is the headspace / liquid volume
ratio r that each rival needs.

**Arithmetic (mine):**
- Air sealed at 25 C and 1 atm holds 40.9 umol gas per mL. At 20.95 % O2 that is **8.56 umol O2 per mL of
  headspace**. Pure O2 holds 40.9 umol/mL.
- Dissolved O2 at sealing is at most 0.26 mmol/L under air and 1.24 mmol/L under O2.
- Every reading charges the precursor at 0.1 mmol/mL. So, per mmol of precursor, the available O2 is
  **85.6 r + 2.6 umol (air)** and **409 r + 12.4 umol (O2-flushed)**.

| rival, arm | O2 needed (umol per mmol) | r needed, air | r needed, O2-flushed |
|---|---|---|---|
| R1, air | 3.8-7.6 | >= 0.06 | any |
| R1, air + Cu | 10-21 | >= 0.21 | n/a (air, fn c) |
| R2, air | 250-700 | **2.9-8.2** | 0.58-1.7 |
| R2, air + Cu | 500-1900 | **5.8-22** | n/a (air, fn c) |

Some example vessels:
- 10 mL in a 20 mL vial: r = 1.0, so 88 umol/mmol under air;
- 1 mL in a 4 mL vial: r = 3, so 259 umol/mmol;
- 1 mL in a 10 mL vial: r = 9, so 773 umol/mmol.

**Outcome: inconclusive on printed data, but sharpened.**
- The O2-flush reading rescues R2's air arm (r >= 0.6). It never rescues the copper arm, because footnote c states
  air.
- At the footnotes' 10 mL scale, a common closed vial (<= 30 mL, r <= 2) refutes R2 under air.
- At the Experimental's 1 mL scale, a vial of >= 7 mL would not refute it.
- So the decisive unknowns are r and the scale.

**Rate flag (assumption-laden, not a verdict).**
- R2's low end consumes 250 umol per 10 mL in 120 min, i.e. 0.21 mmol/L/min.
- Under the sealed air headspace at 100 C (pO2 ~0.26 atm, solubility ~0.75 mmol/L/atm), the liquid holds ~0.2 mmol/L
  O2. So the vessel would need k_L a >= ~1 /min.
- Unstirred closed vials typically run 0.006-0.06 /min.
- P(R2's O2 demand can be transferred in an unstirred closed vial) ~0.2 (mine). R1 needs 0.006 mmol/L/min, which
  any vial supplies.

**Proposed text for MF-33 / M3 open question 11:**

> Headspace unprinted (re-read 2026-11-06). R2 is feasible only if V_headspace / V_liquid >= 2.9 (air arm) and >= 5.8
> (air + Cu arm; footnote c says air). R1 is feasible at any r >= 0.06. The O2-flushed reading of expt II rescues
> R2's air arm (r >= 0.6) but not its copper arm. A design that states r and the scale decides. At the footnotes'
> 10 mL scale in a <= 30 mL vial, R2 is refuted.

## CK-06: the Schibilsky 2019 pots (MF-28)

**Pass rule (inventory, plus mine).** Explained (a pot-definition bug) if the pH is not passed; otherwise the misfit
stands as structural. An integrity error is a unit, MW or concentration conversion in the bundle that does not
reproduce from its quotes.

**Trace (read-only engine calls):**
1. **The pH is passed to the spec.** `panel.core_spec` passes `conditions.ph`, so `spec.process.ph` = 5.0 / 8.0.
   The buffer is "none" (declared unbuffered).
2. **The lane receives it.** Its B15 declaration prints x0.38 (pH 5) / x1.91 (pH 8) on `k_asn_glc` and x0.54 / x1.51
   on `k_acr_dp`. The pot has no asparagine, so neither step carries flux.
3. **The furanics do not respond.** DMHF is 21.84 ug/L and HMF 55409 ug/L at pH 3, 5, 7, 8 and 10 alike.
4. **Routing is the cause.** `engine._PRECURSOR_LANE` sends Ala to the acrylamide lane. The trunk's declared pH law
   is applied only `if lane == TRUNK` (`engine.py`: declarations near l. 1568, parameters near l. 2698). That law
   covers the Amadori-decay steps (Martins 2003, B12) and the B40 3-DG exit term. The acrylamide lane's docstring
   says "No pH dependence".
5. **Lane diagnostic** (the same pot with glycine in place of alanine; a stand-in, not a prediction). In the trunk
   lane:
   - DMHF 1621 / 5158 ug/L (pH 8 / 5 = x3.18; measured 1153 / 5894, x5.1);
   - HMF 20135 / 41390 (x2.06; measured x1.76);
   - the acrylamide lane's DMHF is 74x below the trunk lane's at the same charge.
6. **Bundle arithmetic reproduces.**
   - Charges: 4.504 g / 180.16 / 0.050 L = 0.500 M glucose and 2.228 g / 89.09 / 0.050 L = 0.500 M alanine, diluted
     1:2 to 0.25 M each.
   - Furaneol 9 / 46 uM x 128.13 = 1153.17 / 5893.98 ppb.
   - HMF 454 / 798 uM x 126.11 = 57253.94 / 100635.78 ppb.
   - Furfural 17 / 25 uM x 96.08 = 1633.36 / 2402.0 ppb.

**Outcome:**
- **The pH is passed, and the misfit stands as structural.** Its location is now exact: the acrylamide lane carries
  no pH law on the trunk steps that make HMF and DMHF.
- The source itself could not be re-read, because the Schibilsky thesis is not on the laptop (the bundle says
  "SOURCE NOT ON DISK"). Step 1 (R) for the printed cells stays open until it is downloaded (DepositOnce, TU Berlin
  2019).

**Spent.** The glycine diagnostic was compared with Schibilsky's measured pH ratios. So Schibilsky cannot be the
decisive test of a change that routes alanine pots to the trunk lane or gives the acrylamide lane the trunk's pH law.

**Proposed CF (MF-28).** `kind: engine_structure`, `sources: [schibilsky2019]`:

> The Schibilsky glucose + alanine pots (pH 5 / 8, 130 C, 2 h, unbuffered) pass their pH to the engine. L-alanine
> routes them to the acrylamide lane, whose pH factor acts only on `k_asn_glc` / `k_acr_dp` (no flux without
> asparagine). The trunk's measured Amadori-decay pH law (Martins 2003, B12) and the B40 3-DG exit term are applied
> only in the trunk lane (`engine.py`, `lane == TRUNK`). DMHF and HMF are therefore identical at every pH (21.84 /
> 55409 ug/L). With glycine standing in (trunk lane), the same pot gives DMHF x3.2 and HMF x2.1 from pH 5 to 8
> (measured x5.1 / x1.76) and 74x more DMHF.
>
> Options:
> - (a) run the trunk's declared terms on the trunk steps inside the acrylamide lane;
> - (b) route a non-asparagine amino acid to the trunk lane;
> - (c) declare the acrylamide lane pH-blind for furanics, and refuse.
>
> Shipped outputs move (every acrylamide-lane HMF / DMHF), so the owner decides first. Schibilsky is spent for this
> decision (seen at diagnosis); the decider must be another pH pair.

## CK-02: the zhou2025b xylose-alone rows (MF-21)

**Pass rule (inventory, plus mine).** Explained (integrity / pot definition) if the -16 sigma comes from a zero
prediction or a mis-keyed unit. An integrity error is a keyed value off the plotted marker by more than its read
error, a wrong axis unit, or a wrong charge.

**Evidence:**
- **The four rows are MGO.** In `kinetic_core_joint_deoxyosone_thiol_fit_s042_r3en.json` the data set
  `zhou2025b_fig2a_xylose_alone_dicarbonyls | MGO` (n 4) reads:
  - raw at the optimum, **-0.407 dex**;
  - after its scale, **-1.016 dex**; whitened -16.12;
  - scale "free".

  The same paper's AX-ARP MGO (`zhou2025b_figs1_3_axarp_xylose_120c | MGO`) is **+0.863 dex** raw and +0.254 after
  the scale.
- **The scale is shared.** `series_of` gives Zhou's MGO one free scale per paper x species, so the xylose-alone and
  AX-ARP MGO series share it.
- **The figure agrees with the keyed rows.** In the render of Zhou2025b p. 5:
  - Fig. 2a's axis reads "Concentration (mmol/L)". Its MGO markers sit at ~0.02 / 0.022 / 0.031 / 0.033; keyed
    0.0195 / 0.0216 / 0.0308 / 0.0329.
  - GO at 120 min is keyed 0.1155; the text says "degradation of Xyl itself contributed only 0.12 mmol/L".
  - Fig. 2e's AX-ARP MGO sits at ~0.08 / 0.075 / 0.07 / 0.09; keyed 0.0815 / 0.0758 / 0.0719 / 0.0892.
- **The charge is consistent.** The 20 mM xylose is not printed, but the 120 min dicarbonyls sum to 0.237 mM, which is
  1.18 % of 20 mM against the printed "Only 1.1 %".
- **The conditions agree.** pH 8 set, no buffer, pressure-resistant flask in an oil bath, 120 C, ice-water quench;
  OPD derivatisation, HPLC-UV 315 nm.

**Outcome: integrity holds, and this is not a zero or floor prediction.**
- The -16 sigma is a -1.0 dex residual left after a scale that the AX-ARP pots set, divided by small read errors
  (0.04-0.07 dex).
- The misfit is a **within-paper ratio**: relative to the data, the engine's MGO from AX-ARP is ~1.27 dex too high
  against its MGO from free xylose.
- The laboratory scale (CF-060) cancels inside one paper, so it cannot be the cause.

**Proposed amendment of MF-21's text:**

> Not a zero or floor prediction (raw -0.41 dex, MGO only, 4 rows). The residual is -1.0 dex after the free MGO scale
> shared with the same paper's AX-ARP pots (+0.86 dex raw): a within-paper ratio error of ~1.27 dex between MGO from
> AX-ARP and MGO from free xylose. Integrity holds (re-read 2026-11-06). Candidates: the AX-ARP route to MGO is too
> fast, or free-xylose fragmentation to MGO is too slow. Not the laboratory scale.

## Re-reads of the never-re-read sulfur / Strecker misfits

### MF-10: the Givaudan ribose / NF ratio

**Pass rule (mine).** An integrity error is any of:
- a printed cell that differs from the atlas value;
- a basis that does not cancel in the ratio;
- an analyte other than 2-methyl-3-furanthiol;
- an engine pot that contradicts the printed conditions.

**Evidence:**
- **The table.** de Roos 2005 p. 124, Table 1 (render; the OCR text agrees), "Quantity produced (mg/mole of
  cysteine)". 2-Methylfuran-3-thiol: ribose 43.5, xylose 36, glucose < 0.05, norfuraneol 0.9.
- **The conditions** (footnote a): "Cysteine : carbonyl compound : lactate = 1 : 1 : 1.67; water content: 37% (w/w);
  pH 5.0; boiling for 1.5 h". The data are unpublished (ref. 13, Laarhoven), with no SD, no n and no method.
- **The engine's pot** (route-barrier prereg D2): 100 C, 90 min, pH 5 clamped, ~2 M declared. My concentration check
  agrees: ~1.5 mol/kg of solution, ~1.9 mol/L.
- **Open conditions, already flagged** (`simulable: false`, `condition_uncertain`):
  - the boiling point of a 63 % solids pot, against an SDE pot diluted to ~100 C;
  - the NaOH used for pH 5, which is not counted.

**Outcome: integrity holds.** The ratio is 43.5 / 0.9 = 48.3 per mol of Cys, and the basis cancels.

**Ownership.** CK-40 is owned by M2, and in the parked M1 + M2 prereg v2 Givaudan is a reported guard, not decisive.
This re-read read only cells already SEEN in the atlas (E6) and tested no candidate. List it as SEEN in M2's item 2.

### MF-15: the Grosch pot charges (CK-10)

**Pass rule (CK-10, plus mine).** Explained if the pots carry the wrong charges. An integrity error is a keyed cell
or a reactant amount off the page, or a charge lost between the keyed pot and the engine spec.

**Evidence:**
- **The page agrees with the keyed table** (Grosch & Zeiler-Hilgart 1992 p. 185, Table II and the Experimental):
  - MF-SH: no. 1 4.7 +- 1.2 (N 3), no. 3 5.2 +- 0.4, no. 4 33 +- 3.7 (N 4), no. 6 19 +- 2.5, no. 7 23 +- 1.2,
    no. 8 9.0 +- 2.8 ug.
  - Footnote a: thiamine.HCl 30 umol, TPP 30 umol, ribose 670 umol, cysteine.HCl 325 umol, GSH 230 umol. The masses
    reproduce the moles within 1 %.
  - Conditions: Na-pyrophosphate pH 5.7, 0.2 mol/L, 500 mL; Method I (reflux 2 h, extracted cold); SIDA.
- **The charges reach the engine intact.** `pot_spec` builds thiamine 0.06, cysteine 0.65, ribose 1.34 and
  glutathione 0.46 mmol/L, in 0.2 M pyrophosphate, at 100 C for 120 min, in the sulfur lane.
- **The engine's MFT** (ug/L):

  | pot | no. 1 | no. 3 | no. 4 | no. 6 | no. 7 | no. 8 |
  |---|---|---|---|---|---|---|
  | MFT | 1.0396 | 1.0396 | 1.0398 | 1.1905 | 1.0451 | 1.9201 |

  - Thiamine left: 0.0576 mmol/L in every pot.
  - No. 4 makes H2S (0.0028 mmol/L) and keeps 0.565 mmol/L cysteine.
  - No. 4 / no. 1: measured x7.0, predicted x1.0002.

**Outcome: integrity holds, not explained.** The thiamine route makes MFT first order in thiamine, with no cysteine
or H2S term (candidate (a) of the inventory).

### MF-16: the disulfide share

**Pass rule (mine; CK-09's threshold).** An integrity error is any of:
- the wrong species;
- a basis mismatch between the target and the engine's observable (mass vs molar vs thiol equivalents is 0.30 dex;
  >= 0.5 dex would explain part);
- a keyed value off the page.

**Evidence:**
- **Zhou 2023 Table 1** (Zhou2023b.pdf), bis(2-methyl-3-furyl)disulfide with ARP + Cys: 59.70 / 102.59 / 50.07 ug/L
  at pH 6 / 7 / 8. MFT: 696.99 / 1588.57 / 525.62 ug/L. Authentic standard. All as keyed.
- **Zhang 2024 Fig. 1a** (render). MFT-MFT ("bis (2-methyl-3-furan) disulfides", TCI standard > 98 %):
  - with cysteine, ~0.115 against MFT 1.34 ng/mL;
  - with cystine, 1.09 against 2.01.

  This agrees with the dossier. The temperature is not printed (flagged).
- **The targets match the engine's observable.** The targets (`zhou_pH7_dimer_over_MFT` 0.0323,
  `zhang_fig1_cys_dimer_over_MFT` 0.0429, `zhang_fig1_gcys_dimer_over_MFT` 0.2711) are molar MFTD / MFT, the same
  kind as the engine's ratio. With the true MW of 226.32 they are 0.0326 / 0.0433 / 0.2736: the targets used
  2 x 114.17 = 228.3, an error of +0.004 dex.
- **No 2x slip is possible between the papers.** By-mass and thiol-equivalent shares coincide within 1 %
  (2 x 114.17 / 226.32 = 1.009), so Zhang's by-mass figure and Zhou's thiol-equivalent figure are on the same footing.
- **Measurement lead (step 3, outside this check).** Zhou 2023 SI Table S1 calibrates the disulfide as
  y = 0.0208x - 0.7126, which gives zero response at 34.3 ug/L.
  - With a zero intercept, the same responses read 15.8 / 68.3 / 25.4 ug/L at pH 8 / 7 / 6 (-0.50 / -0.18 /
    -0.37 dex).
  - The true share could therefore be lower, and the engine's x11 deficit at pH 7 would shrink to ~x7. That does not
    close it.
- **Already recorded, not new:** the two share rows are "declared hold-out, in-sample since B60"
  (`fit_m1_supply_combinations.py`).

**Outcome: integrity holds.**

**Proposed correction (immaterial; not under `data/`).** In `generate_kinetic_core_b2_{1,2,3}_fit.py` the targets
should read 0.0323 -> 0.0326, 0.0429 -> 0.0433 and 0.2711 -> 0.2736 (MW of C10H10O2S2 = 226.32).
- Do not edit the frozen generators (the wave manifest pins them).
- Carry the exact values into the next fit's target definitions.
- No shipped output moves: the B60 rows shift 0.004 dex against sigma_log 0.4.

### MF-30: the glucose-pot PA lag (hofmann2000 Fig. 2)

**Pass rule (mine).** An integrity error is any of:
- a point off its marker by more than the read error;
- a wrong axis unit or basis;
- a condition (time zero, temperature, atmosphere) mis-stated in a way that moves the 10-60 min deficit.

**Evidence** (hofmann2000 p. 436, Fig. 2, render):
- The axis reads "amount [umol/mmol glucose]"; the time axis is broken at 250 / 450 min.
- PA sits at ~0.2 / 0.5 / 0.9 / 1.4 / 2.0 / 2.2 and PAA at ~0.1 / 0.3 / 0.8 / 1.8 / 3.2 / 3.5 (10-500 min). The keyed
  values are 0.196 ... 2.204 and 0.096 ... 3.498.
- The recipe (98 C, 1 mmol each in 10 mL, the phosphate print of CF-020) and the duplicate with 2000b Table 2 expt II
  are declared in the keyed file.
- Time zero and the heat-up are not printed.

**Outcome: integrity holds.** A real come-up time would make the measured early PA *lower*. That widens the engine's
-3.04 dex deficit at 10 min, so heat-up cannot explain the lag shape.

### MF-33: the Strecker O2 dependence (source side)

**Pass rule (mine).** An integrity error is a ratio that does not reproduce from the printed cells, or an engine pot
whose identity differs from the printed one without a declared stand-in.

**Evidence:**
- **The ratios reproduce.** PA air / argon is 5.5 / 0.6 = 9.17 (ARP-Phe) and 1.4 / 0.4 = 3.5 (Glc/Phe), an order
  of 2.63x. PA + PAA gives 9.4 for ARP-Phe.
- **B29's engine pots are declared stand-ins** (`kinetic_core_b29_prereg.md` sec. 7): a fed glycine Amadori (DFG) and
  glucose + glycine, with the observable AKG + AKM.
- **The Glc/Phe 3.5 may mix vessels.** It may compare a reflux pot with a closed argon vial (keyed duplicate ruling
  D(i)).
- **Missing deoxygenation can only push the ratio up.** Residual O2 in the argon arm would lower the measured
  air / argon, so the anoxic ratio is >= 9.2. That is further from the engine's ceiling (~1.8-2.8), not closer.

**Outcome: integrity holds.**

---

## What this batch spends

These datasets were read here and become SEEN for any later case. All were already SEEN:
- Hofmann 2000b Table 2, Table 3, the Methods and the Table 2 discussion text;
- Zhou 2025b Fig. 2a / 2e and Secs. 2.3 / 2.5 / 3.2;
- Schibilsky (bundle only), Grosch 1992 Table II, de Roos 2005 Table I, Zhou 2023 Table 1 and SI Table S1, Zhang 2024
  Fig. 1, Hofmann 2000 Fig. 2;
- the Zhang 2025b bundle and its declaration.

Engine outputs were seen for:
- the Schibilsky pots, including the glycine stand-in (so Schibilsky is spent for the MF-28 decision);
- the Grosch Table II pots.

## Owner-level calls this record raises

1. **MF-28.** Choose option (a), (b) or (c). Shipped acrylamide-lane HMF / DMHF would move.
2. **CK-50 / M3.** The R2 verdict needs a stated headspace / liquid ratio and scale. These belong in the design brief
   for DM:strecker_routes (a).
