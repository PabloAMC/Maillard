# M0 integrity batch 1: source re-reads for the trunk, acrylamide, pyrazine and lipid misfits

> **What this is.** Step 1 (integrity) of the misfit protocol
> ([OPERATING_PROCEDURE.md](../../docs/guides/OPERATING_PROCEDURE.md) sec. 4) for eight misfits of the misfit
> inventory (draft PR #82, branch `misfit-inventory` at `12d88121`, `docs/reference/misfit_inventory.md`): MF-31
> (Pan half), MF-39, MF-37, MF-38, MF-42, MF-34, MF-23 and MF-28 (keyed values only). These are checks CK-03, CK-04,
> CK-05 and CK-08 (Pan and Leahy halves) of the inventory's Batch 0. Strecker and sulfur misfits are another agent's batch.
>
> **How.** Each source was re-read from its PDF in the main checkout's article folder (laptop-only; cited here by its
> dossier, never by path), against the keyed rows: conditions, units, basis, analyte, duplicates, digitisation. Figures
> were rendered and read where a scored value sits in one. The corrected residuals are **arithmetic on the recorded
> prediction** (the `core_panel_scores` row or the fit report row). No refit, no forward run, no file under `data/`,
> `src/`, the registers, the inventory, cases or preregs was edited. Built on `origin/main` at `4174ed8d`.
>
> **Blinding.** No blind or reserved dataset was read. Excluded: the Deng 2022 half of MF-31 (M3 phase 2 reserved);
> the E13 / E14 TEST studies and Lee 2024 (MF-34); hofmann2000b Figs. 2-3; B59b hold-outs; WUR. Five of the misfit
> rows (Leitzen "Steinhagen", Lin 2021 "Chang", Ye 2024, Lin 2022, Schibilsky) live in `maillard_path` hold-out bundles.
> Those bundles are scored on the public panel and their values are printed in `core_panel_scores.md`, so they are SEEN.
> Following the inventory's Batch 0 rule, only their methods and the scored cells (plus the control or unheated
> columns of the same tables) were read.
>
> **Order and criteria.** The misfits were re-read largest suspected slip first, in the order given below. The test
> for each one was set before its PDF was opened, and is stated as "error if". It comes from the inventory's
> CK-03/04/05/08 pass rules ("a unit slip of >= 1 dex explains") plus the generic list: units, basis, analyte,
> duplicates, digitisation and pot definition. **Disclosure:** the criteria were held in working notes during the
> reads, and this file was written afterwards. Every outcome is recorded, including "integrity holds".

## Summary

| MF | source (dossier) | outcome | residual as recorded | residual with the correction | what moves it |
|---|---|---|---|---|---|
| MF-31 | Pan 2025 (`pan2025_extraction.md`) | **integrity holds** | -3.6 to -5.6 dex (methional / MeSH rates) | unchanged; the unit-free cross-check at 140 C / 600 s reads **-3.71 dex** | no unit or basis slip can explain it |
| MF-39 | Ma 2024 (`ma2024_extraction.md`) | **integrity error (pot definition)**: the precursor charge is unsourced | -3.98 dex | not computable; the row has no sourced charge. Brackets: -3.86 (formed amount), -3.69 (formed amount and a 37 s residence) | withdraw the charge, then refuse the row |
| MF-37 | Lin 2021 ("Chang 2021"; `lin2021_extraction.md`) | **integrity holds** | +2.34 / +0.40 / +0.65 (AA), -0.44 (HMF) | unchanged | heat-up unprinted: a step-2 lead |
| MF-37 | Ye 2024 (`ye2024_extraction.md`) | **integrity holds** (the 0.2 M molarity stays a labelled inference) | +1.70 dex | +1.70 + log10(c / 0.2 M) (bimolecular, low conversion) | — |
| MF-38 | Lin 2022 (`lin2022_extraction.md`) | **integrity error (prose vs figure)** on HMF; acrylamide holds | AA -0.96, HMF -0.80 | AA -0.96; **HMF -0.41** on the Fig. 4A bar | conflict inside the source |
| MF-42 | Li 2026, Bi 2020, Trikusuma 2020 | **integrity errors (pot definition)**: carrier undiluted (Li, Trikusuma); undeclared carried hexanal (Bi). Values hold | lane median -0.45 | Li **-0.86 to -3.15**; Trikusuma decadienal **+0.21 -> -0.81**; Bi hexanal **-0.47 -> +0.63** | the corrected misfit is deeper and mixed in sign |
| MF-34 | Leahy 1989 (`leahy1989_extraction.md`) | **integrity holds** | -2.88 (B18) / -2.78 (B21) | basis bracket -3.06 / -2.96; lysine's two amines at most +0.30 | — |
| MF-23 | Leitzen 2021 ("Steinhagen"; `leitzen2021_extraction.md`) | **integrity errors (pot definition)**: undeclared unheated levels, end pH keyed as the start pH, vessel unkeyed, F0 not a hold time. Values hold | glucosone -1.80, GO -1.54, HMF -0.97, MGO 0.00 | glucosone **-1.53**, GO **-0.68** (carried GO survives), HMF -0.95, MGO +0.13 | the GO half of the pattern shrinks; glucosone stands |
| MF-28 | Schibilsky 2019 (dissertation) | **not re-read**: the source is not on disk. The keyed conversions check | DMHF -1.72 / -2.43 | unchanged | download needed |

**Net.** Three re-reads hold (MF-31 Pan, MF-34, MF-37) and one is not runnable (MF-28). Four find errors.
- **The errors are in pot definitions, not in keyed values.** Every scored number read matches its print, except
  Lin 2022's HMF, where the print disagrees with itself.
- **Explained by integrity: none.**
  - MF-39 becomes unscorable as keyed.
  - MF-23 loses about 0.9 dex of its glyoxal half.
  - MF-38 loses half of its HMF half.
  - MF-42 gets *deeper* once the carrier dilution is applied (Li, Trikusuma), and its sign becomes mixed (Bi).
- **No misfit here is a unit slip.** The two largest (Pan -3.7, extruded SPI -3.98) are checked against
  unit-independent quantities.

## MF-31 — Pan 2025 B22 rows (Pan half only; Deng reserved)

**Error if:**
- the K unit (unprinted) is mis-inferred by >= 3 dex;
- the methionine charge (mg/L vs mg/mL in the same sentence) or the sugar charge is mis-keyed;
- the keyed rates differ from Table 2;
- the end-point levels are on another basis.

**Re-read** (Methods sec. 2.2; Results sec. 3.1; Table 1; Table 2):
- **Keyed rates.** All nine `PAN_RATES_UMOL_L_MIN` (generator `generate_kinetic_core_b22_fit.py`) are the Table 2
  "K x 10^-4" numbers x 60. They match to the printed digit, e.g. methional 1.823 / 16.8 / 89.9 and MeSH 1.579 / 2.342
  / 6.335.
- **Pot.** "40 mg/L methionine and sugars, including 20 mg/mL fructose, 15 mg/mL glucose, and 15 mg/mL sucrose ... 50
  mM citrate buffer at pH 6.2" is verbatim. The keyed pot is Met 0.268 mM, Fru 111, Glc 83 (83.3), citrate pH 6.2.
  Sucrose is omitted and declared; its ~10 % hydrolysis adds ~2 % to the reducing pool (< 0.01 dex).
- **Methionine in mg/L, not mg/mL.** 40 mg/mL would mean 268 mM Met. The paper's own consumption figures rule that
  out: Met falls "around 8 %" at 120 C / 4 min (21 mM if the charge were 268 mM) while fructose and glucose fall ~4 %
  (~8 mM in all). A 268 mM charge would consume more amine than the reducing sugar lost. At 0.268 mM the paper is
  self-consistent.
- **Unit-independent check.** The text prints "at 140 C for 600 s, methanthiol and methional concentrations peaked at
  26.43 ug/L and 626.31 ug/L".
  - 626.31 ug/L is 6.01 umol/L of methional.
  - B22's report records the engine's pot at 140 C / 10 min as MTAL 1.18e-6 mmol/L (1.18e-3 umol/L).
  - log10(1.18e-3 / 6.01) = **-3.71 dex**. This uses no K, so no unit of K can rescue the row.
  - 626 ug/L sits inside the 1-800 ug/L methional calibration range of Table 1.

**Outcome: integrity holds.** The 3.6-5.6 dex is not a unit or basis slip. The residual stands at -3.71 dex on the
unit-free end point.
- **Side reading (not a test).** The same B22 diagnostic gives the engine's MGO at 140 C / 10 min as 42 uM. Pan prints
  6.2 uM MGO at 120 C / 4 min, a different condition. The donor supply is not obviously short, which points the
  deficit at the methionine Strecker step, candidate (b).
- **Spent:** Pan 2025 Table 2, sec. 3.1 text, Methods (all SEEN).

## MF-39 — extruded soy protein isolate, 130 C (`acrylamide_spi_extrusion_130C_ACSRef3`)

**Error if:**
- the 150 ug/kg is not the Fig. 2D 130 C bar;
- its basis is per dry matter;
- the precursor charge (Asn 42 mM, Glc 42 mM, "SPI 50 mM"), pH 6.0, a_w 0.35 or the 25 s are not in the source.

**Re-read:**
- **Fig. 2D** (rendered, p. 7). The 130 C bar reads ~150 ug/kg (letter a), 150 C ~120, 170 C ~82, control ~38. The
  keyed 150 matches.
- **Methods 3.1-3.2.** The materials are soy protein isolate and corn starch at 9:1. The only chemicals are standards
  and reagents (CML, CEL, acrylamide-d3, MGO, GO, 3-DG, 2,3-BD, OPD). **No asparagine or glucose is added.**
- **Table 1** prints water, protein and total amino acids (TAA, mg/g) only. No free asparagine, no reducing sugar, no pH
  and no water activity appear anywhere.
- **Moisture.** The 130 C product holds 25.0 % water (Table 1) and the feed 30 %. An a_w of 0.35 does not follow from
  either.

**Outcome: integrity error (pot definition).** The measured value holds, but the precursor charge {SPI 50, L-Asn 42,
D-Glc 42 mM}, pH 6.0 and a_w 0.35 are not in Ma 2024. Their origin is the pre-DOI "ACSRef3" intake. Its sister value
(62.62 ppb) was already shown to match nothing in the paper.
- **The -3.98 dex tests a pot the paper does not describe.**
- **Brackets, by arithmetic on the recorded 0.0156 ppb:**
  - formed amount 112 ug/kg (150 - 38 control): **-3.86 dex**;
  - with a 37 s 130 C residence (Unlu 2002 bracket, linear at low conversion): **-3.69 dex**;
  - if the 150 is on a dry basis (unstated), the wet value is 112 ug/kg: -3.86 dex.
- **No charge could rescue the row.** At 42 / 42 mM the engine is 10^3.98 low. The formation is bimolecular at low
  conversion, so closing it would need ~4 M of each precursor, more than the melt holds. Free Asn and reducing sugar in
  an SPI : starch blend are almost certainly far below 42 mM, which only widens the gap.
- **The residual is a conditions problem (step 2), not integrity.** The process is a ten-zone extrusion, keyed as an
  isothermal hold with an unsourced charge. Integrity removes the row; it does not explain it.
- **Spent:** Ma 2024 Fig. 2D, Table 1, sec. 3.1-3.2 (SEEN).

## MF-37 — glucose + Asn, 180 C (Lin 2021 bundles "Chang2021"; Ye 2024)

**Lin 2021. Error if:**
- 28 / 1459 / 832 ppb or HMF 7 ppm are not the AG-arm values at the stated times and vehicles;
- the charge, pH or acetate are mis-keyed;
- the ppb basis is not the reaction solution;
- the 10 and 30 min rows come from different arms.

**Re-read** (sec. 2.1-2.3, 3.1, 3.3; Fig. 4 rendered, p. 7):
- **Acrylamide.** "glucose and asparagine after 10, 20, and 30 min of heating was 28, 912, and 1459 ppb ... (Figure
  3B)" is the AG arm in 1 % acetic acid. The water arm reads 832 ppb at 30 min (Figure 3A,B).
- **HMF.** "7 ppm" is the AG arm. Fig. 4A, AG in water at 30 min, reads ~7.2 ppm. The keyed 7000 ppb matches.
- **Charge.** 0.5 g L-asparagine monohydrate and 0.5 g D-glucose made up to 100 mL, with the pH set to 5.8 and then
  6.0. This gives 33.3 / 27.75 mM, as keyed.
- **Acetate.** 1 % acetic acid is 0.167 M (w/v reading), keyed as 0.17 M.
- **Analysis.** HPLC-UV at 210 nm, calibration 0-3125 ppb; all values lie inside it.
- **Unprinted:**
  - the heating apparatus, vessel, aliquot and come-up time. The sister paper (Lin 2022, same group, same "Chang et
    al." method) names a dry-bath incubator.
  - Lin 2022 cites this paper's Maillard-arm HMF as "5.13 ppm [18]", against the 7 ppm printed here. Fig. 4A shows
    ~5.4 ppm at 20 min and ~7.2 at 30 min. The keyed 30-min value follows the figure.

**Outcome: integrity holds** for every keyed value and condition.
- A unit or basis slip would be shared by the 10 and 30 min rows. It cannot produce a measured x52 rise against an
  engine that falls.
- The candidate is the come-up time (a, step 2). The apparatus is not printed, so it needs a declared band.

**Ye 2024. Error if:**
- 140.58 umol/mol Asn is not the control;
- the denominator is not Asn;
- buffer, vessel or time are mis-keyed.

**Re-read** (sec. 2.2, 3.1, Table 3):
- "The level of acrylamide under the control condition (i.e., 0 g/mol Asn) was 140.58 +/- 13.92 umol/mol Asn". The
  value matches; every dose in the paper is per mol of Asn.
- "Equimolar solution of glucose and asparagine ... phosphate buffer (0.1 M, pH 6.86), 4 mL ... 25 mL thick-walled
  pressurized glass tube in an oil bath at 180 C ... 30 min". This matches.
- The molarity is not printed. The 0.2 M from the cited Knol 2005 method stays a labelled inference.
- Internal-standard arithmetic (1 ug 13C3-AA into 4 mL) is consistent with 0.2 M and also with 0.02 M, so it does not
  decide the molarity.

**Outcome: integrity holds.** The residual stays +1.70 dex at 0.2 M. At low conversion the yield per Asn scales with
[Glc], so the residual moves as +1.70 + log10(c / 0.2 M). It would vanish only at c ~ 4 mM.
- **Spent:** Lin 2021 secs. 2-3, Fig. 3 prose, Fig. 4A/C; Ye 2024 sec. 2.2, sec. 3.1, Table 3 control row (all SEEN).

## MF-38 — fructose + Asn (Lin 2022)

**Error if:**
- 1859 ppb or 12.28 ppm are not the AF-in-water values at 30 min;
- the charge or vehicle is mis-keyed;
- prose and figure disagree.

**Re-read** (sec. 2.2; Figs. 3-4 rendered, pp. 6-7, Fig. 4 at 400 dpi):
- **Acrylamide.** "The heating of asparagine-fructose solution for 30 min generated 1859 and 9401 ppb in water and
  acetic acid, respectively". Fig. 3A shows the AF 30-min bar at ~1860 ppb, with an SD whisker to ~2700 (about +-0.2
  dex). This matches.
- **HMF.** The prose reads "HMF generated by the formation of a dicarbonyl intermediate and 3-deoxyglucosone from the
  Maillard reaction (12.28 ppm (Figure 4A) and 5.13 ppm [18])".
  - The axis was calibrated on the same panel's fructose-only 30-min bar, which reads 91.6 against the printed 91.56.
  - On that axis the **AF 30-min bar in Fig. 4A reads 5.0 ppm** (SD to ~6). In acetic acid (Fig. 4B) it reads ~4.9.
  - No bar in either panel is near 12.28.
  - The same sentence quotes the group's own earlier paper as 5.13 ppm, where that paper prints 7 ppm (Lin 2021 above).
    The prose numbers of this sentence are unreliable.
- **Charge.** Fructose 27.75 mM, Asn monohydrate 33.3 mM, water, pH 6.0, dry-bath incubator. This matches.

**Outcome: integrity error (prose vs figure) on HMF; acrylamide holds.**
- HMF on the figure value: log10(1966 / 5000) = **-0.41 dex**, against -0.80 on the prose value.
- The fructose-axis signature of MF-38 rests on the acrylamide row (-0.96, holds) and half of the HMF row.
- **Spent:** Lin 2022 sec. 2.2, sec. 3.2 prose, Figs. 3A/B and 4A/B (SEEN).

## MF-42 — lipid aldehydes and 2-pentylfuran in matrices (Li 2026, Bi 2020, Trikusuma 2020)

**Error if:**
- a scored value is the wrong row or column;
- the basis (wet / dry, per sample) or the quantification class is mis-keyed;
- an unheated or carried level is omitted;
- the pot's composition, as the engine charges it, differs from the print.

**Re-read:**
- **Li 2026** (sec. 2.3, 2.8, Table 2):
  - **Values.** HMPE-0 min reads hexanal 605.64, nonanal 74.37, 2-pentylfuran 5625.80 ug/kg. This matches.
  - **Basis.** The basis is wet, as far as the print says. The powder is described in sec. 2.5: "HMPEs were ground
    following rapid freezing with liquid nitrogen"; no drying step is printed. Sec. 2.8 uses "2.00 g of powdered HMPE".
    The numbers are semi-quantitative "relative contents" against cyclohexanone, as the bundle declares.
  - **Pot.** "soy protein isolate and wheat gluten hydrolysates at a dry basis ratio of 6:4 ... moisture content
    maintained at approximately 57%", final zone 160 C.
  - The bundle charges "Soy Protein Isolate" alone. The lipid lane turns a carrier into a hydroperoxide pool as
    PV / 2 x lipid fraction x 1 kg/L (`parameters_lipid.looh_charge_mmol_per_l`, reached from
    `engine.py` ~L2410-2437). No mass fraction of the carrier in the food is applied. The engine therefore charges the
    pot as if it were 100 % SPI.
  - The print makes SPI 0.43 x 0.6 = 0.258 of the wet mass. The gluten hydrolysate's lipid is not printed.
  - The bundle's a_w 0.35 contradicts 57 % moisture. The lipid lane reads only temperature and time, so no output moves.
- **Trikusuma 2020** (sec. 2.1-2.2, Table 2):
  - **Values.** UHT 782 / 163 / 24.0 / 46.9 and non-UHT 331 / 59.4 / 8.24 / 0.06 ug/L. This matches, and the carried
    levels are already declared (B31 / B35).
  - **Pot.** "pea protein isolate (PPI, 80% protein, 8% fat, 2% sugar w/w)" at "3% w/w": 0.24 % w/w fat.
  - The engine charges the pea carrier's declared 2.5 % lipid at 100 % of the pot, which is **10.4x (1.02 dex) too
    much fat**.
  - The process includes an 80 C preheat of unprinted duration before the 140 C / 6 s hold, as declared.
- **Bi 2020** (sample preparation, Table 4):
  - **Values.** Roasted peas, hexanal 324 +/- 34.3 ug/kg. This matches.
  - **Raw column.** The same table prints raw peas at **1260 +/- 114** ug/kg hexanal and nonanal at 69.8, with nonanal
    not detected after roasting. The text says "A heat-induced decrease in the content of hexanal".
  - The bundle declares no carried level, so the engine's formed-from-zero 109 ppb is scored against a datum that is
    mostly surviving raw-material hexanal, after a **net 3.9x loss** on roasting.
  - The matrix is whole peas roasted and then milled, proxied by the pea-isolate carrier.

**Outcome: integrity errors (pot definition); every scored value holds.** Corrected residuals, by arithmetic on the
recorded predictions:

| row | recorded | corrected | how |
|---|---|---|---|
| Li 2-pentylfuran | -2.56 | -2.93 to -3.15 | x 0.43 (gluten lipid = SPI lipid) to x 0.258 (gluten lipid 0) |
| Li hexanal | -0.94 | -1.30 to -1.52 | same |
| Li nonanal | -0.49 | -0.86 to -1.08 | same |
| Trikusuma hexanal | -0.34 | -0.37 | formed part x 0.096; carried 331 kept |
| Trikusuma 2-pentylfuran | -0.40 | -0.43 | same |
| Trikusuma nonanal | -0.19 | -0.43 | same |
| Trikusuma (E,E)-2,4-decadienal | +0.21 | **-0.81** | same (almost all formed) |
| Bi hexanal (roasted) | -0.47 | **+0.63** | carried 1260 added, no loss in the lane; matching needs 17 % survival |

- **Integrity does not explain MF-42.** Applying the printed compositions makes the formation rows lower still. The
  one row that flips (Bi) flips because it measures net loss of a carried level, which the lane has no term for. The
  lane-offset statistic (median -0.45, 88 % one sign) is therefore not interpretable until the pots are fixed.
- **Confidence.** P ~ 0.85 that the engine applies no carrier mass fraction. This is from reading the code path; a
  forward run with the lipid fraction scaled would confirm it.
- **Spent:** Li 2026 secs. 2.3, 2.8, Table 2 (HMPE-0 min column); Trikusuma 2020 secs. 2.1-2.2, Table 2; Bi 2020
  sample preparation and Table 4, raw and roasted (all SEEN).

## MF-34 — Leahy 1989 total pyrazines (B18 T4, B21)

**Error if:**
- 13.1 ppm is not the lysine-glucose 2 h / 95 C total;
- the charge, pH or time are mis-keyed;
- ppm is on the 15 mL analysed volume (x1.5);
- the shares do not sum.

**Re-read** (ch. 7 Methods, Table I, Table III, text):
- **Total.** "lysine-glucose system yielded 13.1 ppm total pyrazines while the asparagine-glucose resulted in 0.74
  ppm". This matches `LEAHY_TABLE_III` (55.8 / 41.5 / 2.4 %, total 13.1).
- **Pot.** "0.1M for both amino acid and sugar, in a pH 9.0 0.1M borate buffer"; 10 mL in Teflon-capped tubes. This
  matches the engine pot (Glc 100 + Gly 100 mM, pH 9, 95 C, 120 min). Glycine stands in for lysine, as declared.
- **Internal consistency.** Table I's lysine-glucose pyrazine rate of 3.596 ppm/h x 2 h + 0.06 = 7.25 ppm, against
  Table III's 55.8 % x 13.1 = 7.31 ppm. The two agree.
- **Basis.** It remains open whether ppm refers to the 15 mL sample (10 mL reaction + IS + NaOH) or the reaction. The
  15 mL reading makes the reaction concentration 1.5x higher.

**Outcome: integrity holds.**
- Recorded residuals: -2.88 dex (B18) and -2.78 (B21).
- The basis bracket makes them -3.06 / -2.96.
- Lysine's two amino groups at first order in amine would close at most 0.30 dex.
- The supply deficit stands at step 4 (case M6).
- **Spent:** Leahy 1989 ch. 7 Methods, Tables I and III, Results text (SEEN; not an E13 / E14 TEST study).

## MF-23 — glucose alone, 121 C autoclave (Leitzen 2021, bundle `..._Steinhagen2021`)

**Error if:**
- the six values are not Table 4, scheme A, 121 C;
- the units are not ug/mL of solution;
- 2-KDG is not glucosone;
- the unheated level is omitted;
- pH, time, vessel or O2 are mis-keyed;
- glucose is a monohydrate basis.

**Re-read** (Tables 1, 3, 4, 5, 6; sec. 2.2-2.3; Methods 4.2-4.2.2):
- **Values.** Table 4, scheme A, 121 C: GO 5.6, MGO 2.6, 2-KDG 7.5, 3-DG/3-DGal 52.2, 3,4-DGE 55.5, 5-HMF 17.4 ug/mL,
  n = 27. This matches all six, and 2-KDG is glucosone.
- **Undeclared unheated levels** (Table 6, "not heat sterilized (reference value)"): GO **1.0 +/- 0.5**, MGO **0.9 +/-
  0.1**, 2-KDG **0.1**, 3-DG and 3,4-DGE n.d., 5-HMF 0.1 ug/mL. The bundle charges glucose alone, from zero.
- **pH.** Table 3 gives 4.36 +/- 0.0 for scheme A at 121 C, which is the **end** pH, measured at room temperature
  after cooling (sec. 2.2). The unheated solution is **4.98**.
  - The bundle keys 4.36 as `conditions.ph` with `initial_ph_scale: at_temperature`.
  - Its provenance note says the paper prints "4.4-5.2" and that 4.36 lies below it. That note is wrong: Table 3
    prints 4.36 (scheme A: 5.17 / 4.67 / 4.36 at 111 / 116 / 121 C).
- **Time.** 18 min is F0, the z = 10 C lethality integrated over "heating and cooling phases" (sec. 4.2.2), not a hold
  at 121 C. Scheme B adds the heating and cooling F0 to the hold, so the real cycle is longer than 18 min and its
  profile is unprinted. The bundle runs 18 min isothermal.
- **Vessel.** "225 mL PP bottles (Kautex), each filled with 200 mL glucose solution" (sec. 4.2.1). That is 25 mL of
  air headspace: ~0.21 mmol O2, ~1.1 mM if all of it dissolved. The O2 the measured glucosone (42 uM) and glyoxal
  (97 uM) need is ~0.14 mM. The bundle keys the vessel as "unstated".
- **Glucose.** "glucose monohydrate" was used for a 10 % (w/v) solution per the German monograph. If 10 % refers to
  the monohydrate, the charge is 505 mM, not 555: -0.04 dex.

**Outcome: integrity errors (pot definition); every keyed value holds.** Corrected residuals, by arithmetic on the
panel's recorded predictions, with the unheated levels added and surviving the cycle:

| species | recorded | corrected | note |
|---|---|---|---|
| glucosone | -1.80 | **-1.53** | +100 ppb unheated |
| glyoxal | -1.54 | **-0.68** | +1000 ppb unheated; -1.46 if the carried GO is consumed |
| 5-HMF | -0.97 | -0.95 | +100 ppb |
| methylglyoxal | 0.00 | **+0.13** | +900 ppb |
| 3-DG, 3,4-DGE | +0.10, -0.82 | unchanged | n.d. unheated |

- **Effect on the misfit's signature.** It was "one sign on the oxidative products". After correction the glucosone
  deficit (-1.5) stands, and the glyoxal deficit halves or better.
- **Effects that arithmetic cannot price.** The start pH (4.98, not 4.36) and the longer real cycle both raise every
  product. They need a forward run, which is CK-03's step-2 follow-up and not done here.
- **The O2 arithmetic leaves candidate (a) open.** The vessel held enough O2 for the oxidative products measured.
- **Spent:** Leitzen 2021 Tables 1, 3, 4 (scheme A rows), 5, 6; secs. 2.2-2.3, 4.2-4.2.2 (SEEN). Scheme B and
  Table 5 were seen as context only.

## MF-28 — Schibilsky 2019 furanics (keyed values only; CK-06 is the other agent's)

**Error if:** the uM-to-ppb conversions or the Tabelle 17 rows are mis-keyed.

**Not re-read.** The dissertation (TU Berlin, 2019) is not in the article folder; the bundles say "SOURCE NOT ON DISK;
second-hand". The keyed arithmetic was checked against the quoted Tabelle 17 rows and the thesis's own Tabelle 10
molar masses:

| quantity | quoted | x molar mass | keyed |
|---|---|---|---|
| Furaneol (DMHF) | 9 / 46 uM | 128.13 | 1153.17 / 5893.98 ppb |
| HMF | 454 / 798 uM | 126.11 | 57253.94 / 100635.78 ppb |
| furfural | 17 / 25 uM | 96.08 | 1633.36 / 2402.0 ppb |

All six match. The recorded DMHF residuals reproduce: log10(21.841 / 1153.17) = -1.72 and log10(21.841 / 5893.98) =
-2.43.

**Outcome: integrity not established.** The keyed values are consistent with the quotation, but no primary re-read
was possible. **Download request:** Schibilsky S., "Bildung fluechtiger Verbindungen und Carbonsaeuren in
Modellsystemen unter Maillard-Bedingungen", Dissertation, TU Berlin 2019 (institutional repository). Read Tabelle
17, sec. 5.1.6 / 5.1.6.2 and the vessel choice.
- **Spent:** nothing new.

## Proposed data corrections (NOT applied; for the coordinator's separate PR)

Each correction gives the file, the row, the old value, the new value and the place in the source.

1. **`data/benchmarks/acrylamide_spi_extrusion_130C_ACSRef3.json`** (MF-39).
   - `precursors` {Soy protein isolate 50, L-Asparagine 42, D-Glucose 42 mM}: withdraw. Ma 2024 sec. 3.1 and Table 1
     print no free Asn, no reducing sugar and no molarity.
   - `conditions.ph` 6.0 and `water_activity` 0.35: withdraw (unprinted; Table 1 gives 25 % product water).
   - Then refuse the row ("precursor charge not printed"), per the owner's 2026-11-02 rule for values not in the cited
     paper.
   - **Shipped outputs moved:**
     - `core_panel_scores` (a PRIMARY trust-loop predictive row becomes refused);
     - `core_prediction_uncertainty`, `null_model_baselines`, `data_wishlist`;
     - the JAX acrylamide reference (`jax_trunk_acrylamide_module`);
     - `gap_heatmap`;
     - the frozen B3 fit report (`spi_extrusion_130C_acrylamide_ppb`, sigma 1.0 dex in B3's objective) is not
       regenerated but would be stale, and any B3 refit changes.
   - The shipped acrylamide constants do not move until a refit.
2. **`data/benchmarks/external_validation/maillard_path/mp_holdout_glucose_only_autoclave_121C_Steinhagen2021.json`**
   (MF-23).
   - (a) Add `conditions.carried_volatiles` {glyoxal 1000, methylglyoxal 900, glucosone 100, 5-HMF 100} ppb from
     Leitzen 2021 Table 6. 3-DG and 3,4-DGE are n.d.
   - (b) Change `conditions.ph` from 4.36 (`initial_ph_scale: at_temperature`) to 4.98 (`bench_cooled`, Table 3,
     non-autoclaved). Record 4.36 as the measured end pH (Table 3, scheme A, 121 C), and correct the provenance note
     that claims a printed 4.4-5.2.
   - (c) Change `conditions.vessel` from fill/vessel null and atmosphere "unstated" to fill 200 mL, vessel 225 mL, PP
     bottle (Kautex), air headspace (sec. 4.2.1).
   - (d) Add a note that 18 min is F0 (z = 10 C, heating and cooling included; sec. 4.2.2), not an isothermal hold.
   - **Shipped outputs moved:** `core_panel_scores` rows (glyoxal -1.54 to about -0.68, glucosone -1.80 to -1.53, MGO
     0.00 to +0.13, HMF -0.97 to -0.95, by arithmetic); (b) and (d) need a forward run.
   - Also moved: `core_prediction_uncertainty`, `contrast_scores`, the directional rows that read this pot, and
     `null_model_baselines`.
   - A REFERENCE hold-out: no fit moves.
3. **`data/benchmarks/external_validation/maillard_path/mp_holdout_fructose_asparagine_180C_Lin2022.json`** (MF-38).
   `holdout_targets.5-Hydroxymethylfurfural (HMF).target_value` is 12280 ppb from the prose, which conflicts with the
   Fig. 4A AF 30-min bar (5.0 ppm, SD ~1.0, read at 400 dpi on an axis checked against the printed 91.56).
   - Either replace it with 5000 ppb (`value_provenance: figure_readoff`, uncertainty ~20 %), or keep 12280 with the
     conflict recorded and the row demoted to `ambiguous`.
   - The coordinator decides; proposed CF-P1 below.
   - **Moves:** the `core_panel_scores` HMF row (-0.80 to -0.41), `core_prediction_uncertainty` and `contrast_scores`.
4. **`data/benchmarks/external_validation/external_validation_bi_2020_roasted_pea_hexanal.json`** (MF-42).
   - Add `conditions.carried_volatiles` {hexanal 1260, nonanal 69.8} ppb from Bi 2020 Table 4, raw-pea column, as
     B31 did for Trikusuma.
   - Record that roasting lowered hexanal 3.9x and that nonanal is "—" after roasting.
   - **Moves:** the `core_panel_scores` hexanal row (-0.47 to +0.63 with no loss term), the envelope and the lipid
     lane offset.
   - Better still, refuse the row as "net loss of a carried level; the lane has no loss term" until a volatilisation
     or loss term exists.
5. **The lipid carrier mass fraction** (MF-42; needs a schema field and an engine change, so it is not a pure data
   edit).
   - **`data/benchmarks/pea_isolate_uht_140C_Trikusuma2019.json`:** declare the carrier at 3 % w/w with 8 % fat
     (sec. 2.1-2.2). Today it is implicitly 100 % at 2.5 %.
   - **`data/benchmarks/external_validation/external_validation_li_2026_spi_wg_hme_control.json`:**
     - declare SPI at 0.258 w/w (6:4 dry basis, 57 % moisture; sec. 2.3) plus a wheat-gluten-hydrolysate carrier;
     - withdraw `precursors."Soy Protein Isolate".concentration_mM` 1000 (no molar mass) and `water_activity` 0.35
       (contradicts 57 % moisture).
   - **Moves:** every lipid formation row in `core_panel_scores` (the table in the MF-42 section), the envelope's
     lipid rows, `null_model_baselines`, the lane-offset diagnostic and the W5.4 La cells.
   - This is a CF / choices-register item: see CF-P3.

The Ye 2024 and Lin 2021 bundles, the Pan 2025 generator rows, the Leahy constants and the Schibilsky bundles need no
correction from this batch.

## Proposed CF texts (ids to be assigned by the coordinator; nothing written to the register)

**CF-P1 (conflict, acrylamide lane; sources lin2022, lin2021).**
- **Quantity:** HMF in the asparagine + fructose (and + glucose) Maillard arm, water, pH 6, 180 C, 30 min.
- **Observation:** Lin 2022's prose prints 12.28 ppm "(Figure 4A)" and cites Lin 2021 as 5.13 ppm. Lin 2022 Fig. 4A
  reads 5.0 ppm (fructose-only bar at 91.6 against the printed 91.56 calibrates the axis). Lin 2021 prints 7 ppm, and
  its Fig. 4A reads ~7.2 ppm at 30 min and ~5.4 ppm at 20 min.
- **Explanations:**
  - (i) a prose transcription error in Lin 2022 (12.28 belongs to no bar);
  - (ii) a different time or arm meant.
- **Remedy:** score the figure value, or demote the row; ask the authors.
- **Affects:** `mp_holdout_fructose_asparagine_180C_Lin2022` (HMF).

**CF-P2 (conflict, trunk; source leitzen2021).**
- **Quantity:** the thermal history and starting state of the 121 C autoclave pot.
- **Observation:** the scored 18 min is F0 (z = 10 C lethality including heating and cooling), not a hold; the cycle
  profile is unprinted. The unheated solution already carries GO 1.0, MGO 0.9, glucosone 0.1 and HMF 0.1 ug/mL (Table
  6). The keyed pH 4.36 is the cooled end pH; the start is 4.98. The bottle holds 25 mL of air over 200 mL.
- **Explanations:** MF-23's oxidative deficit is partly a pot-definition artefact (glyoxal about -0.9 dex of its
  -1.54), and the rest is structural (glucosone -1.5 stands).
- **Tests:** forward runs with the start pH and a declared cycle band (CK-03 step 2).
- **Affects:** `mp_holdout_glucose_only_autoclave_121C_Steinhagen2021`.

**CF-P3 (choice / conflict, lipid lane).**
- **Quantity:** the hydroperoxide pool of a matrix pot.
- **Observation:** the lane charges a lipid carrier as the whole pot (PV / 2 x lipid fraction x 1 kg/L), but two
  scored sources print a diluted carrier:
  - Trikusuma: 3 % w/w isolate at 8 % fat, so 10.4x less fat than charged;
  - Li 2026: SPI 0.258 w/w of a 57 %-moisture melt plus an unprinted-lipid gluten hydrolysate.
- **With the printed compositions,** Li's rows fall to -0.9 to -3.2 dex and Trikusuma's decadienal flips from +0.21
  to -0.81.
- **Explanations:** the lane's under-prediction (MF-42) is larger than recorded. Its cause is upstream of the matrix
  layer: the pool's PV, the yields, or carried volatiles in the raw isolates.
- **Remedy:** a carrier mass-fraction field in the bundle schema and engine; re-derive MF-42's size before any
  mechanism work.
- **Affects:** the lipid rows of `core_panel_scores`.

**CF-P4 (conflict, lipid lane; source bi2020).**
- **Quantity:** hexanal on roasting peas, 160 C / 30 min.
- **Observation:** raw 1260 +/- 114, roasted 324 +/- 34.3 ug/kg (Table 4), a net 3.9x loss. The panel scores the
  roasted value as formation from zero.
- **Explanations:** volatilisation in the open oven, and consumption of the carried hexanal.
- **Remedy:** declare the carried level, and refuse or score as loss once the lane has a loss term.
- **Affects:** `external_validation_bi_2020_roasted_pea_hexanal`.

## Process note

Ten curated JSON files (the bundles above and the intake registry) were printed for reading with the system `python3`,
not through `./scripts/docker_maillard.sh run`, which AGENTS.md requires. These were read-only prints; no number in
this record was computed by them. All arithmetic above was run in the `maillard` env (native runtime). PDF text was
extracted with `pdftotext` and figures rendered with `pdftoppm` into the session scratchpad; nothing was written under
`data/`.
